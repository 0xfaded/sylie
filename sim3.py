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
import se3
import common

def v(w):
    l = w[6]
    w = w[:3,:]

    if len(l.free_symbols) == 0 and abs(l.evalf()) < 1e-9:
        return se3.v(w[:3,:])

    s = skew(w)
    e = sy.exp(l)

    theta = norm3(w)
    l2t2 = l**2 + theta**2

    A = (e-1)/l

    if len(theta.free_symbols) == 0 and abs(theta.evalf()) < 1e-9:
        B = (1 + e*(l-1)) / l**2
        C = (-2 + e*(2 - 2*l + l**2))/(2*l**3)
    else:
        B = (theta*(1 - e*sy.cos(theta)) + e*l*sy.sin(theta)) / (theta*l2t2)
        C = A/(theta**2) - e*sy.sin(theta)/(theta*l2t2) - l*(e*sy.cos(theta)-1)/(theta**2*l2t2)

    return A*sy.eye(3) + B*s + C*s*s

def dv(w, dw_dx):
    l = w[6]
    w = w[:3,:]

    if len(l.free_symbols) == 0 and abs(l.evalf()) < 1e-9:
        return se3.dv(w, dw_dx)

    dw_dx = dw_dx[:3,:]

    wW = (w.T * dw_dx)[0,0]
    
    s = skew(w)
    S = skew(dw_dx)
    e = sy.exp(l)


    theta = norm3(w)
    l2t2 = l**2 + theta**2

    e = sy.exp(l)
    le1cos = l * (1 - e * sy.cos(theta))

    F1 = -2*(l*e*sy.sin(theta)) / (theta*l2t2**2) - \
        2*(1 - e*sy.cos(theta)) / (l2t2**2) + \
        (l*e*sy.cos(theta) - e*sy.cos(theta) + 1) / (theta**2 * l2t2) + \
        (e*sy.sin(theta)) / (theta * l2t2) - \
        l*e*sy.sin(theta) / (theta**3 * l2t2) - \
        (1 - e*sy.cos(theta)) / (theta**2*l2t2)

    F2 = -2*le1cos / (theta ** 2 * l2t2 ** 2) - \
        2*le1cos / (theta**4 * l2t2) + \
        2*e*sy.sin(theta) / (theta * l2t2**2) - \
        e*sy.cos(theta) / (theta**2 * l2t2) + \
        l*e*sy.sin(theta) / (theta**3 * l2t2) + \
        e*sy.sin(theta) / (theta**3 * l2t2) - \
        2*(e-1) / (l*theta**4)

    F3 = l*e*sy.sin(theta) / (theta * l2t2) + \
        (1 - e*sy.cos(theta)) / l2t2

    F4 = le1cos / (theta**2 * l2t2) - \
        e * sy.sin(theta) / (theta * l2t2) - \
        (-e + 1) / (l * theta**2)

    return wW * (F1 * s + F2 * s*s) + F3 * S + F4 * (s*S + S*s)

def dv_dlambda(w):
    l = w[6]

    w = w[:3,:]
    s = skew(w)
    e = sy.exp(l)

    theta = norm3(w)
    l2t2 = l**2 + theta**2

    if len(l.free_symbols) == 0 and abs(l.evalf()) < 1e-9:
        A = 1/2

        if len(theta.free_symbols) == 0 and theta.evalf() < 1e-9:
            B = 1/3
            C = 1/8
        else:
            B = -(theta*sy.cos(theta) - sy.sin(theta))/theta**3
            C = -(-theta**2 + 2*theta*sy.sin(theta) + 2*sy.cos(theta) - 2)/(2*theta**4)

    else:
        A = e/l - (e-1)/l**2

        if len(theta.free_symbols) == 0 and theta.evalf() < 1e-9:
            B = (l**2*e - 2*l*e + 2*e - 2)/l**3
            C = (l**3*e - 3*l**2*e + 6*l*e - 6*e + 6)/(2*l**4)

        else:
            B = -2*l*(l*e*sy.sin(theta) + theta*(-e*sy.cos(theta) + 1))/(theta*(l2t2)**2) + \
                (l*e*sy.sin(theta) - theta*e*sy.cos(theta) + e*sy.sin(theta))/(theta*(l2t2))

            C = 2*l**2*(e*sy.cos(theta) - 1)/(theta**2*(l2t2)**2) + \
                2*l*e*sy.sin(theta)/(theta*(l2t2)**2) - \
                l*e*sy.cos(theta)/(theta**2*(l2t2)) - \
                e*sy.sin(theta)/(theta*(l2t2)) - \
                (e*sy.cos(theta) - 1)/(theta**2*(l2t2)) + \
                e/(l*theta**2) - \
                (e - 1)/(l**2*theta**2)

    return A*sy.eye(3) + B*s + C*s*s

def exp(w):
    result = sy.eye(4)

    result[:3,:3] = so3.exp(w)
    result[3,3] = sy.exp(-w[6])
    result[:3,3] = v(w) * w[3:6,:]

    return result

def log(M):
    result = sy.zeros(7,1)

    result[:3,:] = so3.log(M[:3,:3])
    result[6] = -sy.ln(M[3,3])
    result[3:6,:] = v(result) **-1 * M[:3,3]

    return result

def dw_dw(w):
    result = sy.zeros(7,7)

    vw = v(w)
    vi = vw ** -1
    so3_dw_dw = so3.dw_dw(w[:3,:])

    result[:3,:3] = so3_dw_dw
    result[3:6,3:6] = vi * sy.exp(-w[6])
    result[6,6] = 1

    for i in range(3):
        dx = sy.zeros(3)
        dx[i] = 1
        
        result[3:6,i] = -vi * (dv(w, so3_dw_dw[:,i]) - skew(dx) * vw) * w[3:6,:]

    result[3:6,6] = -vi * dv_dlambda(w) * w[3:6,:]

    return result

def dba_da(wa, wb):
    return -dba_db(wa, wb)

def dba_db(wa, wb):
    A = exp(wa) ** -1
    B = exp(wb)
    AB = A * B
    w = log(AB)

    vb = v(wb)
    vba_i = v(w) ** -1

    db_db = dw_dw(wb)
    dba_db = so3.dba_db(wa[:3,:], wb[:3,:])

    result = sy.zeros(7,7)
    result[:3,:3] = dba_db
    result[3:6,3:6] = vba_i*A[:3,:3]*B[3,3]
    result[6,6] = 1

    for i in range(3):
        dvba_db = dv(w, dba_db[:,i])
        dvb_db = dv(wb, db_db[:3,i])

        result[3:6, i] = vba_i * (-dvba_db * w[3:6,:] + so3.exp(-wa[:3,:]) * (dvb_db * wb[3:6,:] + vb * db_db[3:6,i]))

    dvba_dlambda = dv_dlambda(w)
    dvb_dlambda = dv_dlambda(wb)

    result[3:6, 6] = vba_i * (-dvba_dlambda * w[3:6,:] + \
            so3.exp(-wa[:3,:]) * (dvb_dlambda * wb[3:6,:] + vb * db_db[3:6,6]) - B[3,3] * A[:3, 3])

    return result
