
import sympy

C_SYMBOL = sympy.symbols("c", commutative=True)
P_SYMBOL = sympy.symbols("p", commutative=True, integer=True)
N_SYMBOL = sympy.symbols("n", commutative=True, integer=True)
Q_SYMBOL = sympy.symbols("q", commutative=True, integer=True)
MU_SYMBOL = sympy.symbols("mu", commutative=True)
TAU_SYMBOL = sympy.symbols("tau", commutative=True)
THETA_SYMBOL = sympy.symbols("theta", commutative=True)
ALPHA_SYMBOL = sympy.symbols("α", commutative=True)

K_SYMBOL_FOR_INPUT_VALIDATION = sympy.symbols("k", commutative=True)


# -1 means infinite on localhost.
REQUEST_TIMEOUT__LONG__LOCALHOST = 10000
REQUEST_TIMEOUT__SHORT__LOCALHOST = 10000

REQUEST_TIMEOUT__LONG__PRODUCTION = 60.0
REQUEST_TIMEOUT__SHORT__PRODUCTION = 5.0

def get_request_timeout_in_secs(path, is_localhost):
    request_timeout = -1
        
    if path == "/validate_phi_psi_input":
        if is_localhost == True:
            request_timeout = REQUEST_TIMEOUT__SHORT__LOCALHOST
        else:
            request_timeout = REQUEST_TIMEOUT__SHORT__PRODUCTION
    elif path.startswith("/module_cruncher"):
        if is_localhost == True:
            request_timeout = REQUEST_TIMEOUT__LONG__LOCALHOST
        else:
            request_timeout = REQUEST_TIMEOUT__LONG__PRODUCTION
    
    return request_timeout


DISALLOWED_SYMBOLS_IN_INPUT = [
    C_SYMBOL, Q_SYMBOL, MU_SYMBOL, TAU_SYMBOL, K_SYMBOL_FOR_INPUT_VALIDATION
]

TAU_SYMBOL_NON_COMM = sympy.symbols("tau", commutative=False)
C_SYMBOL_NON_COMM = sympy.symbols("c", commutative=False)
MU_SYMBOL_NON_COMM = sympy.symbols("mu", commutative=False)

A_SYMBOL = sympy.symbols("a", commutative=True)
A_SYMBOL_NON_COMM = sympy.symbols("a", commutative=False)
B_SYMBOL_FOR_TESTS = sympy.symbols("b", commutative=True)

K_EXPONENT = sympy.symbols("k", commutative=True, integer=True, positive=True)


COEFFICIENT_OF_ONE = 1
TAU_DEGREE_OF_ZERO = 0;

# Set to -1 for normal use. May want to set to a positive number if hitting an infinite recursion case, and the first, or second, etc. reduction step reproduces the issue
# and you want to stop there.
# Infinite recursion most likely means that sympy has to be told a little more explicitly what to do, e.g. some
# simplification of a certain symbolic expression has to have its hand held a little more.
FEATURE_FLAG__REDUCTION_LIMIT = -1
FEATURE_FLAG__REDUCTION_LIMIT_OTHER = -1
# FEATURE_FLAG__REDUCTION_LIMIT_EXT = 10
# FEATURE_FLAG__REDUCTION_LIMIT_OTHER = 1

FEATURE_FLAG__USE_REQUEST_CACHE = True
FEATURE_FLAG__RUN_UNIT_TESTS = False
FEATURE_FLAG__STOP_AFTER_UNIT_TESTS = True
FEATURE_FLAG__INCLUDE_INTERMEDIATE_REDUCTION_MATRICES = False

# NOTE: This is only partially implemented so can't just set to true.
# For all integers it is assumed that e.g. `10^q*tau = 10*tau = tau*10``.
FEATURE_FLAG__APPLY_TWIST_TO_NUMBERS = False


# Pretty much settled on True now, though it makes client response parsing a little more complicated.
# On the flip side it allows this python server to drop abandoned connections and not take up CPU resources running
# calculations if e.g. the client closed the browser.
FEATURE_FLAG__USE_JSON_CHUNKING = True

# Only relevant if FEATURE_FLAG__USE_JSON_CHUNKING == True
MINIMUM_SECONDS_BETWEEN_CHUNKED_CONNECTION_LIVENESS_CHECKS = 5;
