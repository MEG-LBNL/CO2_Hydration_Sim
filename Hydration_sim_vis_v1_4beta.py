#Output intermediate concentrations for the system at steady-state.
from Hydration_sim_main_v1 import *
import matplotlib.pyplot as mt
import matplotlib.ticker as ticker
from scipy.interpolate import griddata
mt.rcParams.update({'font.size': 14, 'lines.markersize': 16})

#Simple procedure comparing kinetics between reversible and classic Michaelis-Menten enzyme rate models.
def compare_kinetics(CO2aq, pH): #CO2aq concentration variable is a vector of values.
    Ka1 = 4.45e-7; Ka2 = 4.69e-11; DIC = 1.8e-3; H = 10**-pH
    H2CO3 = DIC*H**2/(H**2 + H*Ka1 + Ka1*Ka2)
    k2a = HCAII(CO2aq, H2CO3)[2]; k2b = MM(CO2aq)
    print H2CO3, '[H2CO3]  / M'
    figure = mt.gcf(); figure.set_size_inches(8,8)
    mt.subplot(2, 1, 1)
    mt.plot(CO2aq, k2a*Et, 'r+', label = 'revised CA enzyme model')
    mt.plot(CO2aq, k2b*Et, 'c-', linewidth = 2.0, label = 'classical MM model')
    #mt.axhline(y = kcat*Et, linestyle = '--', color = 'k')#maximum enzyme velocity (MM model)
    mt.ylabel('rate, M sec$^-$$^1$'); mt.legend()
    mt.subplot(2, 1, 2)
    residue = abs((k2a - k2b)*Et)
    mt.plot(CO2aq, residue, 'm-')
    mt.xlabel('[ CO$_{2(aq)}$ ]'); mt.ylabel('|Residual|, M sec$^-$$^1$')
    mt.tight_layout()
    mt.show()

#Calculate CO2 hydration and mineralization contours.
def simulate(pH_low, pH_high): #Relevant only where dproton = 0
    endstates = []; enzymerates  = []
    CO2 = np.arange(1.0e-5, 2.1e-5, 1e-6)
    for each in tqdm(CO2):
        for every in np.arange(pH_low, pH_high, 0.1):#Range of aqueous CO2 concentrations.
            H = 10**-every
            H2CO3 = DICi*H**2/(H**2 + H*Ka1 + Ka1*Ka2)
            HCO3 = DICi*H*Ka1/(H**2 + H*Ka1 + Ka1*Ka2)
            CO3 = DICi*Ka1*Ka2/(H**2 + H*Ka1 + Ka1*Ka2)
            init_state[7] = H; init_state[2] = each
            states = odeint(ODE_system, init_state, t, args=(kcat, km))
            endstates.append(states[-1])
            kcat1 = HCAII(each, H2CO3)[2]
            enzymerates.append(kcat1*Et)
    endstates = np.array(endstates); enzymerates = np.array(enzymerates)
    #state matrix has units of concentration; need to multiply by reactor volume processed each sec. to get mol quantity
    flux_volume = 0.01 #seawater flux in a 1 second interval; 6.94e-7 for a 60 ml volume in 24 hrs.
    endstates[:,3] = endstates[:,3]*flux_volume #CaCO3 accumulation by end of simulation time.
    endstates[:,8] = endstates[:,8]*flux_volume #Ca(OH)2 accumulation by end of simulation time.

    print('\n Rendering plot...\n')
    print(min(enzymerates), max(enzymerates))
    mt.style.use('bmh')
    # Plot total sequestered carbonate as f(n) of total [CO2aq], solution pH
    figure = mt.gcf(); figure.set_size_inches(14, 8.0)
    
    # Formatter for two decimal places
    formatter = ticker.FuncFormatter(lambda x, pos: '{:.1e}'.format(x))
    decimal_format = ticker.FuncFormatter(lambda x, pos: '{:.2f}'.format(x))
    
    # Plot total sequestered carbonate as f(n) of total [CO2aq], solution pH
    mt.subplot(2, 3, 1)
    z = endstates[:, 3]
    levels = np.linspace(min(z), max(z), 300)
    x = 1e6 * (endstates[:, 2]); y = -np.log10(endstates[:, 7])
    X, Y = np.meshgrid(x, y); Z = griddata((x, y), z, (X, Y), method='nearest')
    cf = mt.contourf(X, Y, Z, levels, cmap='jet')
    mt.plot(15, 9, 'k.') #plot experimental data
    cbar_ticks = np.linspace(min(z), max(z), 5)
    cbar = mt.colorbar(cf, label='Mol CaCO$_3$', format=formatter, ticks=cbar_ticks)
    #cbar.ax.yaxis.set_major_formatter(formatter)
    mt.ylabel('pH')
    
    # Plot total amount of metal hydroxide formed as f(n) of total [CO2aq], solution pH
    mt.subplot(2, 3, 2)
    z = endstates[:, 8]
    levels2 = np.linspace(min(z), max(z), 300)
    x = 1e6 * endstates[:, 2]; y = -np.log10(endstates[:, 7])
    X, Y = np.meshgrid(x, y); Z = griddata((x, y), z, (X, Y), method='nearest')
    cf = mt.contourf(X, Y, Z, levels2, cmap='jet')
    mt.plot(15, 9, 'k.') #plot experimental data
    cbar_ticks = np.linspace(min(z), max(z), 5)
    cbar = mt.colorbar(cf, label='Mol Ca(OH)$_2$', format=formatter, ticks = cbar_ticks)
    #cbar.ax.yaxis.set_major_formatter(formatter)
    # mt.xlabel('[CO$_{2 (aq)}$] / micromolar')
    # mt.ylabel('pH')
    
    # Plot enzyme catalysis rate as f(n) of pH and aqueous [CO2]
    mt.subplot(2, 3, 3)
    z = enzymerates*1e6 #units in micromolar
    levels3 = np.linspace(min(z), max(z), 300)
    x = 1e6 * endstates[:, 2]; y = -np.log10(endstates[:, 7])
    X, Y = np.meshgrid(x, y); Z = griddata((x, y), z, (X, Y), method='nearest')
    cf = mt.contourf(X, Y, Z, levels3, cmap='jet')
    mt.plot(15, 9, 'k.') #plot experimental data
    cbar_ticks = np.linspace(min(z), max(z), 5)
    cbar = mt.colorbar(cf, label='Catalyzed Reaction Rate, $\mu$M s$^-$$^1$', ticks = cbar_ticks)
    #cbar.ax.yaxis.set_major_formatter(formatter)
    # mt.xlabel('[CO$_{2 (aq)}$] / micromolar')
    # mt.ylabel('pH')
    
    # Plot conversion to H2CO3
    mt.subplot(2, 3, 4)
    z = endstates[:, 4] * flux_volume
    levels = np.linspace(min(z), max(z), 300)
    x = 1e6 * (endstates[:, 2]); y = -np.log10(endstates[:, 7])
    X, Y = np.meshgrid(x, y); Z = griddata((x, y), z, (X, Y), method='nearest')
    cf = mt.contourf(X, Y, Z, levels, cmap='jet')
    mt.plot(15, 9, 'k.') #plot experimental data
    cbar_ticks = np.linspace(min(z), max(z), 5)
    cbar = mt.colorbar(cf, label='Mol H$_2$CO$_3$', format=formatter, ticks=cbar_ticks)
    #cbar.ax.yaxis.set_major_formatter(formatter)
    mt.xlabel('[CO$_{2 (aq)}$] / $\mu$M')
    mt.ylabel('pH')
    
    # Plot conversion to HCO3-
    mt.subplot(2, 3, 5)
    z = endstates[:, 5] * flux_volume
    levels = np.linspace(min(z), max(z), 300)
    x = 1e6 * (endstates[:, 2]); y = -np.log10(endstates[:, 7])
    X, Y = np.meshgrid(x, y); Z = griddata((x, y), z, (X, Y), method='nearest')
    cf = mt.contourf(X, Y, Z, levels, cmap='jet')
    mt.plot(15, 9, 'k.') #plot experimental data
    cbar_ticks = np.linspace(min(z), max(z), 5)
    cbar = mt.colorbar(cf, label='Mol HCO$_3$$^-$', format=formatter, ticks=cbar_ticks)
    #cbar.ax.yaxis.set_major_formatter(formatter)
    mt.xlabel('[CO$_{2 (aq)}$] / $\mu$M')
    # mt.ylabel('pH')
    
    # Plot conversion to CO32-
    mt.subplot(2, 3, 6)
    z = endstates[:, 6] * flux_volume
    levels = np.linspace(min(z), max(z), 300)
    x = 1e6 * (endstates[:, 2]); y = -np.log10(endstates[:, 7])
    X, Y = np.meshgrid(x, y); Z = griddata((x, y), z, (X, Y), method='nearest')
    cf = mt.contourf(X, Y, Z, levels, cmap='jet')
    mt.plot(15, 9, 'k.') #plot experimental data
    cbar_ticks = np.linspace(min(z), max(z), 5)
    cbar = mt.colorbar(cf, label='Mol CO$_3$$^2$$^-$', format=formatter, ticks=cbar_ticks)
    #cbar.ax.yaxis.set_major_formatter(formatter)
    mt.xlabel('[CO$_{2 (aq)}$] / $\mu$M')
    # mt.ylabel('pH')
    
    mt.tight_layout()
    mt.show()

#Generate a contour plot of carbonate conversion as a fn of enzyme kcat and Km; kcat and Kms are vector inputs.
kcats = np.arange(0, 5e6, 5e5); Kms = np.arange(5e-3, 1.8e-2, 1e-3)
#enzyme_optimize(kcats, Kms, 9.0, 1.5e-5)
def enzyme_optimize(kcats, Kms, pH, CO2): #Relevant only where dproton = 0
    endstates = [] #run calculation for Case 1
    CA_params = []
    H = 10**-pH; init_state[7] = H; init_state[2] = CO2
    for kcat in kcats:
        for km in Kms:#Range of aqueous CO2 concentrations.
            states = odeint(ODE_system, init_state, t, args=(kcat, km))
            endstates.append(states[-1]); CA_params.append([kcat, km])
            H2CO3 = DICi*H**2/(H**2 + H*Ka1 + Ka1*Ka2); kcat1 = HCAII(CO2, H)[2]
    endstates = np.array(endstates); CA_params = np.array(CA_params)
    
    #state matrix has units of concentration; need to multiply by reactor volume processed each sec. to get mol quantity
    flux_volume = 0.01
    endstates[:,3] = endstates[:,3]*flux_volume*1e3
    
    # Formatter for two decimal places
    formatter = ticker.FuncFormatter(lambda x, pos: '{:.2e}'.format(x))
    decimal_format = ticker.FuncFormatter(lambda x, pos: '{:.0f}'.format(x))
    
    #Plot total sequestered carbonate as f(n) of total [CO2aq], solution pH
    figure = mt.gcf(); figure.set_size_inches(7.5, 10) #old = (14.4, 5.4)
    mt.style.use('bmh')
    mt.subplot(2,1,1)
    levels = np.linspace(min(endstates[:,3]), max(endstates[:,3]), 1000)
    x = CA_params[:,0]/1e6; y = CA_params[:,1]*1000; z = endstates[:,3]
    X,Y = np.meshgrid(x,y); Z = griddata((x,y), z, (X,Y), method = 'nearest')
    cf = mt.contourf(X, Y, Z, levels, cmap='jet')
    cbar_ticks = np.linspace(min(z), max(z), 5)
    cbar = mt.colorbar(cf, label='mMol CaCO$_3$ Hr$^-$$^1$', format=decimal_format, ticks = cbar_ticks)
    mt.plot(4.4, 12.5, 'c.', label = 'S. Azorense CA')
    mt.plot(1.4, 9.3, 'k.', label = 'HCA II') #
    mt.plot(0.38, 15.8, 'w.', label = 'HCA IX') #https://doi.org/10.1006/bbrc.2001.5824; at pH 8.5, 25 C.
    #mt.colorbar(label = 'Mol Carbonate day$^-$$^1$')#, ticks = [0, 1e7, 2e7, 3e7, 4e7, 5e7, 6e7]) #, ticks = [0, 1e7, 2e7, 3e7, 4e7, 5e7, 6e7]
    mt.xlabel('Enzyme k$_{cat}$ / 10$^6$ [s$^-$$^1$]'); mt.ylabel('Enzyme K$_m$ [mM]')
    mt.tight_layout()
    mt.legend()
    
    mt.subplot(2,1,2)
    efficiency = CA_params[:,0]/CA_params[:,1] #kcat/Km
    mt.plot(efficiency / 1e9, endstates[:,3], 'k+')
    mt.xlabel('Enzyme efficiency, k$_{cat}$ / K$_{m}$ [M$^-$$^1$ s$^-$$^1$] / 10$^9$')
    mt.ylabel('mMol CaCO$_3$ Hr$^-$$^1$')
    mt.tight_layout(); mt.show()