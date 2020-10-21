# no cut, cut
def diff(a, b):
	return 100 * abs(a - b) / b

print("Sensitivity 1, percentage comparison between both cuts and no cut")
print("TOTAL FLUX", diff(1.27508891111732e-09, 1.794113925439598e-09))

print("Wavelet 1")
print("Ngce", diff(5716286.446456, 8500868.807031))
print("Nr", diff(16.381114, 30.549597))
print("Rr", diff(0.092844, 0.113540))

print("Wavelet 2")
print("Ngce", diff(147153.023271, 218835.874846))
print("Nr", diff(57.020094, 97.635438))
print("Rr", diff(0.330510, 0.379332))

print("GLC")
print("Ngce", diff(447.876122, 666.050624))
print("Nr", diff(78.735255, 124.390420))
print("Rr", diff(0.687501, 0.724644))

print("GCE")
print("Ngce", diff(23258.429108, 34588.339187))
print("Nr", diff(7.183457, 19.631184))
print("Rr", diff(0.033646, 0.059368))

print("AIC")
print("Ngce", diff(244980.433781, 364318.084366))
print("Nr", diff(4.713956, 11.860771))
print("Rr", diff(0.025058, 0.039459))

print("NPTF")
print("Ngce", diff(683.051931, 965.281423))
print("Nr", diff(25.805357, 111.243629))
print("Rr", diff(0.066103, 0.255508))

print("Disk")
print("Ngce", diff(17477.628840, 25991.529854))
print("Nr", diff(14.130952, 30.215623))
print("Rr", diff(0.106138, 0.134195))


print("\nSensitivity 1, percentage comparison between resolved cut and no cut")

print("Wavelet 1")
print("Ngce", diff(8500868.807031, 8500868.807031))
print("Nr", diff(11.943565, 30.549597))
print("Rr", diff(0.048110, 0.113540))

print("Wavelet 2")
print("Ngce", diff(218835.874846, 218835.874846))
print("Nr", diff(41.573682, 97.635438))
print("Rr", diff(0.171264, 0.379332))

print("GLC")
print("Ngce", diff(666.050624, 666.050624))
print("Nr", diff(57.406331, 124.390420))
print("Rr", diff(0.356250, 0.724644))

print("GCE")
print("Ngce", diff(34588.339187, 34588.339187))
print("Nr", diff(5.237500, 19.631184))
print("Rr", diff(0.017435, 0.059368))

print("AIC")
print("Ngce", diff(364318.084366, 364318.084366))
print("Nr", diff(3.436973, 11.860771))
print("Rr", diff(0.012984, 0.039459))

print("NPTF")
print("Ngce", diff(965.281423, 965.281423))
print("Nr", diff(23.771308, 111.243629))
print("Rr", diff(0.044567, 0.255508))

print("Disk")
print("Ngce", diff(25991.529854, 25991.529854))
print("Nr", diff(10.302959, 30.215623))
print("Rr", diff(0.054999, 0.134195))
