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

from common import skew, norm3
import so3
import common

def v(w):
    w = w[:3,:]

    s = skew(w)
    theta = norm3(w)

    if len(theta.free_symbols) == 0 and theta.evalf() < 1e-9:
        A = 1/2
        B = 1/6
    else:
        A = (1 - sy.cos(theta)) / theta ** 2
        B = (theta - sy.sin(theta)) / theta ** 3

    return sy.eye(3) + A*s + B*s*s

def dv(w, dw_dx):
    w = w[:3,:]
    dw_dx = dw_dx[:3,:]

    wW = (w.T * dw_dx)[0,0]
    
    s = skew(w)
    S = skew(dw_dx)

    theta = norm3(w)

    F1 = sy.sin(theta) / theta ** 3 - 2 * (1-sy.cos(theta)) / theta ** 4
    F2 = (1 - sy.cos(theta)) / theta ** 4 - 3 * (theta - sy.sin(theta)) / theta ** 5
    F3 = (1 - sy.cos(theta)) / theta ** 2
    F4 = (theta - sy.sin(theta)) / theta ** 3

    return wW * (F1 * s + F2 * s*s) + F3 * S + F4 * (s*S + S*s)


def exp(w):
    result = sy.eye(4)

    result[:3,:3] = so3.exp(w)
    result[:3,3] = v(w[:3,:]) * w[3:6,:]

    return result

def log(M):
    result = sy.zeros(6,1)

    result[:3,:] = so3.log(M[:3,:3])
    result[3:,:] = v(result[:3,:]) **-1 * M[:3,3]

    return result

def dw_dw(w):
    result = sy.zeros(6,6)

    vw = v(w)
    vi = so3.dw_dw(w[:3,:])

    result[:3,:3] = result[3:,3:] = vi

    for i in range(3):
        dx = sy.zeros(3)
        dx[i] = 1
        
        result[3:,i] = -vi * (dv(w[:3,:], vi[:,i]) - skew(dx) * vw) * w[3:6,:]

    return result

def dba_da(wa, wb):
    return -dba_db(wa, wb)

def dba_db(wa, wb):
    AB = exp(wa) ** -1 * exp(wb)
    w = log(AB)

    vb = v(wb)
    vba_i = v(w) ** -1

    db_db = dw_dw(wb)
    dba_db = so3.dba_db(wa[:3,:], wb[:3,:])

    result = sy.zeros(6,6)
    result[:3,:3] = result[3:,3:] = dba_db

    for i in range(3):
        dvba_db = dv(w, dba_db[:,i])
        dvb_db = dv(wb, db_db[:3,i])

        result[3:, i] = vba_i * (-dvba_db * w[3:,:] + so3.exp(-wa[:3,:]) * (dvb_db * wb[3:,:] + vb * db_db[3:,i]))

    return result
