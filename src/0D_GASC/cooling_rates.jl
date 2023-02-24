def Lint_ADAS(Tmin, Tmax,Texp,Lexp, Zimp; N = 100):
# Perform weighted cooling rate integral over specified temperature interval
# Inputs:  Tmin   minimum temperature for integral (eV)
#          Tmax   maximum temperature for integral (eV)
#          Texp   Exponent for temperature weighting
#          Lexp   Exponent for cooling rate weighting
#          Zimp   Z of impurity (Ar: 18; Kr: 36; Xe: 54)

   
    adas_data = load_adas_cooling_rates(Zimp)
        
    dT = LinRange(Tmin,Tmax,N)
    Lz = np.interp(Te,adasdata[Zimp]['temp'],adasdata[Zimp]['Lztot'])
    return np.trapz((Te^Texp)*(Lz^Lexp),Te) *1.e-6
end