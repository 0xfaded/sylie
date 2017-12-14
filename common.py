# BSD 3-Clause License
# 
# Copyright (c) 2017, Carl Chatfield
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# 
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# 
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


import sympy as sy

def skew(w):
    return sy.Matrix([
        [0, -w[2], w[1]],
        [w[2], 0, -w[0]],
        [-w[1], w[0], 0]])


def norm3(w):
    return sy.sqrt(w[0]**2 + w[1]**2 + w[2]**2)

def dodgy_subs(x, theta):
    inv_trig = [
            (sy.acos(sy.cos(theta)), theta),
            (sy.asin(sy.sin(theta)), theta)]

    return x.replace(sy.Abs, sy.Id).subs(inv_trig)

def validate_dw_dw(dw_dw, exp, log, dim, w=None):
    if not w:
        w = sy.randMatrix(dim, 1, -1e6, 1e6) * 1e-6

    numeric = sy.zeros(dim, dim)

    for i in range(dim):
        dx = sy.zeros(dim, 1)
        dx[i] = 1e-6

        w2 = log(exp(dx) * exp(w))

        numeric[:,i] = (w2 - w) / 1e-6

    symbolic = dw_dw(w)

    check_matrices(numeric, symbolic, w, None)

def validate_dba_da(dba_da, exp, log, dim, a=None, b=None):
    if not a:
        a = sy.randMatrix(dim, 1, -1e6, 1e6) * 1e-6
    if not b:
        b = sy.randMatrix(dim, 1, -1e6, 1e6) * 1e-6

    numeric = sy.zeros(dim, dim)

    for i in range(dim):
        dx = sy.zeros(dim, 1)
        dx[i] = 1e-6

        w1 = log(exp(a)**-1 * exp(b))
        w2 = log((exp(dx) * exp(a))**-1 * exp(b))

        numeric[:,i] = (w2 - w1) / 1e-6

    symbolic = dba_da(a, b)

    check_matrices(numeric, symbolic, a, b)

def validate_dba_db(dba_db, exp, log, dim, a=None, b=None):
    if not a:
        a = sy.randMatrix(dim, 1, -1e6, 1e6) * 1e-6
    if not b:
        b = sy.randMatrix(dim, 1, -1e6, 1e6) * 1e-6

    numeric = sy.zeros(dim, dim)

    for i in range(dim):
        dx = sy.zeros(dim, 1)
        dx[i] = 1e-6

        w1 = log(exp(a)**-1 * exp(b))
        w2 = log(exp(a)**-1 * exp(dx) * exp(b))

        numeric[:,i] = (w2 - w1) / 1e-6

    symbolic = dba_db(a, b)

    check_matrices(numeric, symbolic, a, b)

def check_matrices(numeric, symbolic, a, b):
    dim = numeric.cols
    for i in range(dim):
        for j in range(dim):
            if abs(numeric[i,j] - symbolic[i,j]) > 1e-4:
                print('a')
                sy.pprint(a)
                if b:
                    print('\nb')
                    sy.pprint(b)
                print('\nnumeric')
                sy.pprint(numeric)
                print('\nsymbolic')
                sy.pprint(symbolic)
                print()
                raise Exception('numeric[{0},{1}] != symbolic[{0},{1}]'.format(i,j))

