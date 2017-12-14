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

import so3, se3, sim3
import common

for i in range(10):
    common.validate_dw_dw(so3.dw_dw, so3.exp, so3.log, 3)
    common.validate_dba_db(so3.dba_db, so3.exp, so3.log, 3)
    common.validate_dba_da(so3.dba_da, so3.exp, so3.log, 3)

    common.validate_dw_dw(se3.dw_dw, se3.exp, se3.log, 6)
    common.validate_dba_db(se3.dba_db, se3.exp, se3.log, 6)
    common.validate_dba_da(se3.dba_da, se3.exp, se3.log, 6)

    common.validate_dw_dw(sim3.dw_dw, sim3.exp, sim3.log, 7)
    common.validate_dba_db(sim3.dba_db, sim3.exp, sim3.log, 7)
    common.validate_dba_da(sim3.dba_da, sim3.exp, sim3.log, 7)
