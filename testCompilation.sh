#!/bin/bash
# Tests compilation etc for different use cases

# Basic run then Himalaya run then high accuracy mode
autoreconf -vif && ./configure CPPFLAGS="-I/usr/include/eigen3" --enable-two-loop-gauge-yukawa-compilation --enable-himalaya && make -j4 && ./softpoint.x leshouches < inOutFiles/lesHouchesInput && ./configure --enable-two-loop-gauge-yukawa-compilation --enable-three-loop-rge-compilation && make -j4 && ./softpoint.x sugra --tol=1.0e-5 --m0=1000 --m12=800 --a0=-1000 --tanBeta=50 --sgnMu=1 --mt=173.2 --alpha_s=0.1187 --mbmb=4.18 --two-loop-gauge-yukawa --three-loop-rges && ./configure && make -j4
