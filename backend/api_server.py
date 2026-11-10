import sys
sys.dont_write_bytecode = True

from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import logging
import json
import math

import sympy
from sympy import Pow, Mul, Symbol, Poly, Add, Matrix, MatMul, symbols, zeros, ones, simplify, latex, expand, degree, Wild, Number
from functools import reduce

import traceback
import tee_view
import tee_tests
import tee_utils
import tee_phi_psi
import tee_ai_bridge
import tee_consts

import time

from tee_consts import *
from tee_types import *

REQUEST_CACHE__PHI_PSI_EMBEDDING = ThreadSafeLRUCache(1024)

# Just in case some sympy.simplify() calls go a bit haywire.
# UPDATE: Figured out the source of the issue I was dealing with.
# See usage of tee_utils.simplify_each_term_in_expression
# Setting the recursion limit high just seems to make it take longer
# to throw the inevitable recursion limit exception.
# sys.setrecursionlimit(5000)


def parse_phi_psi_options(request_json):
    show_reduction_steps = request_json.get("show_reduction_steps", False)
    replace_c_with_tau = request_json.get("replace_c_with_tau", True)
    replace_q_with_integer = request_json.get("replace_q_with_integer", True)
    group_coefficients = request_json.get("group_coefficients", True)
    strip_integer_exponents_for_symbolic_execution = request_json.get("strip_integer_exponents_for_symbolic_execution", True)
    execute_symbolically = request_json.get("execute_symbolically", True)
    prime_characteristic = int(request_json.get("prime_characteristic", -1))
    extension_degree = int(request_json.get("extension_degree", -1))
    reduce_log_output = False
    only_generate_ext = request_json.get("only_generate_ext", False)

    phi_psi_options = PhiPsiOptions(show_reduction_steps, replace_c_with_tau, replace_q_with_integer, group_coefficients, strip_integer_exponents_for_symbolic_execution, execute_symbolically, prime_characteristic, extension_degree, reduce_log_output, only_generate_ext)

    return phi_psi_options

def handle_phi_psi_input_validation_request(request_json):
    Phi_latex = request_json.get("phi_module")
    Psi_latex = request_json.get("psi_module")

    phi_psi_options = parse_phi_psi_options(request_json)

    phi_psi_input = tee_phi_psi.validate_input(Phi_latex, Psi_latex, phi_psi_options)

    response_dict = tee_view.phi_psi_input_to_dict(phi_psi_input)

    return response_dict

def handle_phi_psi_embedding_request(request_json):
    unit_test_results = None
    if FEATURE_FLAG__RUN_UNIT_TESTS:
        unit_test_results = tee_tests.run_tests()

    if FEATURE_FLAG__RUN_UNIT_TESTS and FEATURE_FLAG__STOP_AFTER_UNIT_TESTS:
        response_dict = {
            "error": "Running unit tests only, turn off FEATURE_FLAG__STOP_AFTER_UNIT_TESTS if needed."
        }
        
        return response_dict

    phi_psi_options = parse_phi_psi_options(request_json)

    start_time = int(time.time_ns())
    Phi_latex = request_json.get("phi_module")
    Psi_latex = request_json.get("psi_module")
    phi_psi_input = tee_phi_psi.validate_input(Phi_latex, Psi_latex, phi_psi_options)

    input_validation_end_time = int(time.time_ns())
    input_validation_elapsed_time = (input_validation_end_time-start_time) / 1e9


    ai_reduction: PhiPsiReductionOutputFromAiCode = tee_ai_bridge.run_reduction(phi_psi_input, phi_psi_options)
    # print(f"ai_reduction={ai_reduction}")

    # phi_psi_reduction_output: PhiPsiReductionOutput = tee_phi_psi.run_analysis(phi_psi_input, phi_psi_options)

    # print(f"STARTING MATRIX COMPARISON")
    # tee_ai_bridge.compare_matrices("W_Matrix", phi_psi_reduction_output.Ext_matrix.raw, ai_reduction.Ext_Matrix.raw, phi_psi_options)
    # tee_ai_bridge.compare_matrices("V_Matrix", phi_psi_reduction_output.V_matrix.raw, ai_reduction.V_matrix)
    # tee_ai_bridge.compare_matrices("W_V_Matrix", phi_psi_reduction_output.W_V_matrix.raw, ai_reduction.W_V_matrix)
    # print(f"ENDING MATRIX COMPARISON")
    

    # response_dict = phi_psi_reduction_output.response_dict
    response_dict = ai_reduction.response_dict

    analysis_end_time = int(time.time_ns())
    analysis_elapsed_time = (analysis_end_time-start_time) / 1e9

    response_dict["unit_test_results"] = unit_test_results
    response_dict["seconds_to_validate_input"] = input_validation_elapsed_time
    response_dict["total_seconds_to_generate"] = analysis_elapsed_time
    
    return response_dict


class SimpleAPIHandler(BaseHTTPRequestHandler):

    # NOTE: Timeouts here only apply to request parsing.
    # Timing out the overall request needs additinal logic elsehwere as well.
    # Realistically since we're not uploading huge files or anything, should never
    # get timeouts triggered here.
    def parse_request(self):
        if not super().parse_request():
            return False  # Parsing failed, stop execution

        is_localhost = tee_utils.is_running_on_localhost(self)
        request_timeout = tee_consts.get_request_timeout_in_secs(self.path, is_localhost)
        
        if request_timeout > 0:
            self.connection.settimeout(request_timeout)
            
        return True

    def do_OPTIONS(self):
        self.send_response(200)

        # Allow all origins
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'GET, POST, PUT, DELETE, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', 'Content-Type, Authorization')
        self.send_header('Access-Control-Max-Age', '86400')

        self.end_headers()

    def do_GET(self):

        self.send_response(200)
        self.send_header('Content-type', 'text/html')
        self.end_headers()

        self.wfile.write("GET request for {}".encode('utf-8'))

        print(f"GET not allowed")
    def do_POST(self):

        is_localhost = tee_utils.is_running_on_localhost(self)
        print(f"is_running_on_localhost={is_localhost}")

        current_epoch_seconds = int(time.time())
        tee_utils.thread_local.response = self
        tee_utils.thread_local.is_localhost = is_localhost
        tee_utils.thread_local.processing_start_time = current_epoch_seconds
        
        if self.path == "/validate_phi_psi_input":
            print(f"Handling {self.path} request...")

            catch_all_exception_message = ""

            content_length = int(self.headers['Content-Length'])
            post_data = self.rfile.read(content_length) 
            post_data_string = post_data.decode('utf-8')

            request_json = json.loads(post_data_string)

            self.send_response(200)
            self.send_header('Access-Control-Allow-Origin', '*')
            self.send_header('Access-Control-Allow-Methods', 'GET, POST, PUT, DELETE, OPTIONS')
            self.send_header('Access-Control-Allow-Headers', 'Content-Type, Authorization')
            self.send_header('Access-Control-Max-Age', '86400')
            self.send_header('Content-type', 'application/json')
            self.end_headers()
            
            try:
                response_dict = handle_phi_psi_input_validation_request(request_json)
                response_dict["success"] = True

            except Exception as e:

                traceback.print_exc()
                catch_all_exception_message = f"Error: {e}"

                response_dict = {
                    "error": catch_all_exception_message
                }
            
            response_bytes_to_return = tee_utils.get_final_json_chunk_bytes(response_dict)
            
            tee_utils.write_final_json_chunk(self, response_bytes_to_return, use_json_chunking=False)
            
        elif self.path == "/module_cruncher":

            print(f"Handling request...")

            catch_all_exception_message = ""

            content_length = int(self.headers['Content-Length'])
            post_data = self.rfile.read(content_length) 
            post_data_string = post_data.decode('utf-8')

            request_json = json.loads(post_data_string)

            self.send_response(200)
            self.send_header('Access-Control-Allow-Origin', '*')
            self.send_header('Access-Control-Allow-Methods', 'GET, POST, PUT, DELETE, OPTIONS')
            self.send_header('Access-Control-Allow-Headers', 'Content-Type, Authorization')
            self.send_header('Access-Control-Max-Age', '86400')
            self.send_header('Content-type', 'application/json')
            # self.send_header("Content-Length", len(response_bytes))
            if FEATURE_FLAG__USE_JSON_CHUNKING:
                self.send_header("Transfer-Encoding", "chunked")
            self.end_headers()

            variant_raw = request_json.get("variant")

            request_cache = None
            request_cache_key = ""

            if FEATURE_FLAG__USE_REQUEST_CACHE == True:
                if variant_raw == Variant.VARIANT__PHI_PSI_EMBEDDING.name:
                    request_cache = REQUEST_CACHE__PHI_PSI_EMBEDDING
                elif variant_raw == Variant.VARIANT__RIGHTMOST_BASIS.name:
                    request_cache = REQUEST_CACHE__RIGHTMOST_BASIS

            response_bytes = None

            if request_cache is not None:
                request_cache_key = tee_utils.unique_hash(post_data_string)
                response_bytes = request_cache.get(request_cache_key)

            if response_bytes is not None:
                print("Returning cached result.")
                tee_utils.write_final_json_chunk(self, response_bytes)
            else:
                try:
                    if variant_raw == Variant.VARIANT__PHI_PSI_EMBEDDING.name:
                        response_dict = handle_phi_psi_embedding_request(request_json)
                    else:
                        raise ValueError(f"Unhandled variant: {variant_raw}")

                    if "error" not in response_dict:
                        response_dict["success"] = True

                except Exception as e:

                    traceback.print_exc()
                    catch_all_exception_message = f"Error: {e}"

                    response_dict = {
                        "error": catch_all_exception_message
                    }
                
                response_bytes_to_return = tee_utils.get_final_json_chunk_bytes(response_dict)
                
                response_dict["cached"] = True
                response_bytes_to_cache = tee_utils.get_final_json_chunk_bytes(response_dict)

                if request_cache is not None and response_dict.get("success") == True:
                    request_cache.put(request_cache_key, response_bytes_to_cache)
                
                tee_utils.write_final_json_chunk(self, response_bytes_to_return)
        else:
            print(f"Invalid path")
            self.send_response(404)

def run(server_class=ThreadingHTTPServer, handler_class=SimpleAPIHandler, port=8080):
    logging.basicConfig(level=logging.INFO)
    server_address = ('0.0.0.0', port)
    httpd = server_class(server_address, handler_class)
    
    logging.info('Starting httpd on port %s...\n', port)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    httpd.server_close()
    logging.info('Stopping httpd...\n')

if __name__ == '__main__':
    from sys import argv

    if len(argv) == 2:
        run(port=int(argv[1]))
    else:
        run()
