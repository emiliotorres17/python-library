



# U is the field that the derivative is being taken of.
# dim is the number of grid points in each direction in U.




# Setup Wavenumber Space
k = np.fft.fftfreq(dim) * dim
Kfield = np.array(np.meshgrid(k,k,k,indexing='ij'))

# LES Level
St[0] = 0.5*np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(U[0]) + 1j*Kfield[0]*np.fft.fftn(U[0])).real
St[1] = 0.5*np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(U[1]) + 1j*Kfield[1]*np.fft.fftn(U[0])).real
St[2] = 0.5*np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(U[2]) + 1j*Kfield[2]*np.fft.fftn(U[0])).real
St[3] = 0.5*np.fft.ifftn(1j*Kfield[1]*np.fft.fftn(U[1]) + 1j*Kfield[1]*np.fft.fftn(U[1])).real
St[4] = 0.5*np.fft.ifftn(1j*Kfield[1]*np.fft.fftn(U[2]) + 1j*Kfield[2]*np.fft.fftn(U[1])).real
St[5] = 0.5*np.fft.ifftn(1j*Kfield[2]*np.fft.fftn(U[2]) + 1j*Kfield[2]*np.fft.fftn(U[2])).real
Rt[0] = 0.5*np.fft.ifftn(1j*Kfield[1]*np.fft.fftn(U[0]) - 1j*Kfield[0]*np.fft.fftn(U[1])).real
Rt[1] = 0.5*np.fft.ifftn(1j*Kfield[2]*np.fft.fftn(U[1]) - 1j*Kfield[1]*np.fft.fftn(U[2])).real
Rt[2] = 0.5*np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(U[2]) - 1j*Kfield[2]*np.fft.fftn(U[0])).real
