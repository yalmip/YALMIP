function top = startofSDPCone(K)

top = 1 + K.f + K.l + sum(K.q) + 3*K.e + sum(K.p);