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

from common import skew, norm3, dodgy_subs
import common

def exp(w):
    s = skew(w)
    theta = norm3(w)

    if len(theta.free_symbols) == 0 and theta.evalf() < 1e-9:
        A = 1
        B = 1/2
    else:
        A = sy.sin(theta) / theta
        B = (1 - sy.cos(theta)) / theta ** 2

    return sy.eye(3) + A*s + B*s*s

def log(R):
    RRT = (R - R.T) / 2
    theta = sy.acos((R.trace() - 1) / 2)

    if len(theta.free_symbols) == 0 and theta.evalf() < 1e-9:
        logR = RRT
    else:
        logR = theta / sy.sin(theta) * RRT

    alpha = logR[2,1]
    beta = logR[0,2]
    gamma = logR[1,0]

    return sy.Matrix([[alpha],[beta],[gamma]])

def dw_dw(w):
    s = skew(w)
    theta = norm3(w)

    A = -1/2
    B = 1/theta**2 - (1 + sy.cos(theta))/(2*theta*sy.sin(theta))

    return sy.eye(3) + A*s + B*s*s

def dba_da(wa, wb):
    return -dba_db(wa, wb)

def dba_db(wa, wb):
    A = exp(-wa)
    B = exp(wb)

    AB = A * B

    theta = sy.acos((AB.trace() - 1) / 2)

    s1 = (sy.sin(theta) - theta * sy.cos(theta)) / (2*sy.sin(theta) ** 2) * (AB - AB.T)

    result = sy.zeros(3,3)

    for i in range(3):
        dx = sy.zeros(3,1)
        dx[i] = 1.0

        dAXB = A * skew(dx) * B
        dtheta = -dAXB.trace() / (2*sy.sin(theta))

        foo = A * exp(1e-6 * dx) * B
        theta2 = sy.acos((foo.trace() - 1) / 2)

        s2 = theta / (2*sy.sin(theta)) * (dAXB - dAXB.T)

        dalpha = dtheta * s1[2,1] + s2[2,1]
        dbeta = dtheta * s1[0,2] + s2[0,2]
        dgamma = dtheta * s1[1,0] + s2[1,0]

        result[0,i] = dalpha
        result[1,i] = dbeta
        result[2,i] = dgamma

    return result
