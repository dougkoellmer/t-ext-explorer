"""Two complete examples. Run with ordinary Python, without SageMath."""
import sympy as sp
from fj_sympy import (
    Parameters, ReductionData, prime_field, function_field, FJExactSequence)


def main():
    # Step 1: fix the exact coefficient field and q-Frobenius.
    theta, tau = sp.symbols("θ τ")
    field = function_field(prime_field(3), theta)
    parameters = Parameters(field, q=3, theta=theta, tau=tau)

    # Step 2: enter the original presentations; Phi is never conjugated.
    examples = {
        "running": (
            [[theta+tau**2, 0], [tau**3, theta+tau**2+tau**3]],
            [[theta+tau]]),
        "fourteen-dimensional": (
            [[theta+tau**4, tau], [0, theta+tau**3]],
            [[theta+theta**2*tau, 0], [1+tau, theta+theta*tau]])}

    for name, (phi, psi) in examples.items():
        # Step 3: search for a jump basis and validate the reduction pivots.
        data = ReductionData(phi, psi, parameters)
        # Step 4: compute both actions, inclusion, projection and certificates.
        report = FJExactSequence(data)
        print(name)
        print("sigma:", report["Sigma"])
        print("jump degrees:", report["JumpDegrees"])
        print("dimensions (M, M0, h):",
              (report["FullDimension"], report["ExtTau0Dimension"], report["H"]))
        print("exact sequence verified:", report["ExactSequenceVerified"])
        for key in ("Pi", "Pi0", "T", "P"):
            print(key, "=", report[key].display())


if __name__ == "__main__":
    main()
