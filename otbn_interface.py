import struct
from typing import List, Tuple
from .polynomials import PolynomialRing
from hw.ip.otbn.util.otbn_sim_py import run_sim

POLY_BYTES = 256 * 4


def dmem_to_coeffs(raw_dmem: bytes, offset: int) -> List[int]:
    coeffs = []
    for w in struct.iter_unpack("<32s", raw_dmem[offset:offset + POLY_BYTES]):
        # parse as signed 32-bit integer
        for v in struct.iter_unpack("<i", w[0]):
            assert len(v) == 1
            coeffs.append(v[0])
    return coeffs


def to_ntt_otbn(poly: PolynomialRing.Polynomial) -> PolynomialRing.Polynomial:
    coeffs_bytes = struct.pack(f"{len(poly.coeffs)}i", *poly.coeffs)
    _, raw_dmem = run_sim("ntt_dilithium_test", [(0, coeffs_bytes)])

    poly.coeffs = dmem_to_coeffs(raw_dmem, 1024)
    return poly


def from_ntt_otbn(poly: PolynomialRing.Polynomial) -> PolynomialRing.Polynomial:
    coeffs_bytes = struct.pack(f"{len(poly.coeffs)}i", *poly.coeffs)
    _, raw_dmem = run_sim("intt_dilithium_test", [(0, coeffs_bytes)])

    poly.coeffs = dmem_to_coeffs(raw_dmem, 0)
    return poly


def ntt_coefficient_multiplication_otbn(poly1: List[int], poly2: List[int], q: int) -> PolynomialRing.Polynomial:
    coeffs1_bytes = struct.pack(f"{len(poly1)}i", *poly1)
    coeffs2_bytes = struct.pack(f"{len(poly2)}i", *poly2)
    _, raw_dmem = run_sim("poly_pointwise_dilithium_test", [(0, coeffs1_bytes), (1024, coeffs2_bytes)])

    new_coeffs = dmem_to_coeffs(raw_dmem, 2048)
    return new_coeffs


def poly_add_otbn(poly1_coeffs: List[int], poly2_coeffs: List[int]) -> List[int]:
    coeffs1_bytes = struct.pack(f"{len(poly1_coeffs)}i", *poly1_coeffs)
    coeffs2_bytes = struct.pack(f"{len(poly2_coeffs)}i", *poly2_coeffs)
    _, raw_dmem = run_sim("poly_add_dilithium_test", [(0, coeffs1_bytes), (1024, coeffs2_bytes)])

    new_coeffs = dmem_to_coeffs(raw_dmem, 2048)
    return new_coeffs


def poly_sub_otbn(poly1_coeffs: List[int], poly2_coeffs: List[int]) -> List[int]:
    coeffs1_bytes = struct.pack(f"{len(poly1_coeffs)}i", *poly1_coeffs)
    coeffs2_bytes = struct.pack(f"{len(poly2_coeffs)}i", *poly2_coeffs)
    _, raw_dmem = run_sim("poly_sub_dilithium_test", [(0, coeffs1_bytes), (1024, coeffs2_bytes)])

    new_coeffs = dmem_to_coeffs(raw_dmem, 2048)
    return new_coeffs


def poly_reduce32_otbn(poly1_coeffs: List[int]) -> List[int]:
    coeffs1_bytes = struct.pack(f"{len(poly1_coeffs)}i", *poly1_coeffs)
    _, raw_dmem = run_sim("poly_reduce32_dilithium_test", [(0, coeffs1_bytes)])

    new_coeffs = dmem_to_coeffs(raw_dmem, 0)
    return new_coeffs


def key_pair_otbn(zeta: bytes) -> Tuple[bytes, bytes]:
    # skip the first 8 bytes (=model, op_state)
    _, raw_dmem = run_sim("key_pair_dilithium_test", [(43840, zeta)]) #43840

    pk = raw_dmem[40000:40000 + 1312]
    sk = raw_dmem[40000 + 1312:40000 + 1312 + 2528]

    return pk, sk


def verify_otbn(pk_bytes: bytes, m: bytes, sig_bytes: bytes) -> int:
    # skip the first 8 bytes (=model, op_state)
    # account for alignment of message
    pk_addr = 40000
    sig_addr = pk_addr + 1312
    m_addr = (sig_addr + 2416)
    m_addr = (m_addr + (m_addr % 32))  # align
    m_len_addr = m_addr + 9 * 4 + 3300
    m_len_bytes = len(m).to_bytes(4, "little")
    print(f"m_len_addr {m_len_addr}")
    print(f"m_len_bytes {m_len_bytes.hex()}")
    regs, _ = run_sim("verify_dilithium_test", [(pk_addr, pk_bytes),
                                                (sig_addr, sig_bytes),
                                                (m_addr, m),
                                                (m_len_addr, m_len_bytes)])

    print(regs)

    # a0 is 0 on success, -1 on fail
    return regs["x10"] == 0

def sign_otbn(sk_bytes: bytes, m: bytes) -> int:
    # skip the first 8 bytes (=model, op_state)
    # account for alignment of message
    sk_addr = 52800
    sig_addr = sk_addr + 2528
    m_addr = sig_addr + 2420
    from math import ceil
    m_addr = int(ceil(m_addr / 32) * 32)  # align
    m_len_addr = m_addr + 9 * 4 + 3300
    m_len_bytes = len(m).to_bytes(4, "little")
    print(f"m_len_addr {m_len_addr}")
    print(f"m_len_bytes {m_len_bytes.hex()}")
    regs, raw_dmem = run_sim("sign_dilithium_test", [(sk_addr, sk_bytes),
                                                (m_addr, m),
                                                (m_len_addr, m_len_bytes)])
    sig = raw_dmem[sig_addr:sig_addr+2420]

    print(regs)

    # sig, 0, siglen
    return sig, regs["x10"], regs["x11"]
