
print("P_tot = ", P_tot)
print("boundary_heat_flux_uniform = ", boundary_heat_flux_uniform)
print("boundary_heat_flux_nonuniform = ", boundary_heat_flux_nonuniform)
print("relative error uniform = ",100*abs(P_tot - boundary_heat_flux_uniform)/P_tot)
print("relative error nonuniform = ",100*abs(P_tot - boundary_heat_flux_nonuniform)/P_tot)