import matplotlib
matplotlib.rc("text", usetex=True)
import matplotlib.pyplot as plt
import matplotlib.animation as animation
plt.switch_backend("agg")
import numpy as np

nfiles = 100
first = 0
interval = 1
fileprefix   = 'checkpoint_'
filesuffix   = '.dat'
saveasprefix = 'sn'

G        = 6.67e-8
pc       = 3.0857e18 # cm
yr       = 31557600.0 # sec
Rsun     = 7.0e10
Msun     = 2.0e33
cmtokm   = 1.0e5
vmax     = 2.0e9
rhoISM   = 1.6e-24
pressISM = 1.0e-5*rhoISM*vmax*vmax
T_Start  = 0.0
T_End    = 40.0e10
vComp    = 1.28e8 # cm/s
mComp    = 1.0*Msun
radComp  = 0.1*Rsun
t0       = 25.0 # yr
gamma    = 5.0/3.0
Runiv    = 8.314e7 # cgs
xmax     = 15.0
xminmc   = 1.0e-4
xmaxmc   = 1.0e3

T = T_End - T_Start
dt = T/float(nfiles)

time = np.zeros(nfiles)
filename = []
for i in range(0, nfiles) :
    num = 10000 + first + i * interval
    numstr = str(num)
    cut = numstr[1:7]
    filename.append(cut)
    time[i] = dt*float(first+i) + t0*yr

rhoComp   = np.zeros(nfiles)
pressComp = np.zeros(nfiles)
velComp   = np.zeros(nfiles)
XComp     = np.zeros(nfiles)
TComp     = np.zeros(nfiles)
csComp    = np.zeros(nfiles)
j = 0
for name in filename :

    print 'starting',name
    data = np.loadtxt(fileprefix + name + filesuffix)

    r     = data[:,0]
    dr    = data[:,1]
    rho   = data[:,2]
    press = data[:,3]
    v     = data[:,4]
    X     = data[:,5]

    menc = 0.0
    mcoord = np.zeros(len(rho))
    for k in range(0,len(rho)) :
        mshell = rho[k]*4.0*np.pi*r[k]*r[k]*dr[k]
        menc = menc + mshell
        mcoord[k] = menc

    rComp = vComp*time[j]/pc
    cs2   = gamma*press/rho
    cs    = np.sqrt(cs2)
    T     = cs2/gamma/Runiv
    #mach = np.sqrt(v*v/cs2)

    rin   = r < rComp*pc
    m_in  = rho[rin]*4.0*np.pi*r[rin]*r[rin]*dr[rin]
    mcoordComp = m_in.sum()/Msun

    compIndex = rin.sum() # upper
    upperIndex = compIndex+1
    lowerIndex = compIndex-1
    rhoComp[j] = rho[lowerIndex:upperIndex].mean()
    pressComp[j] = press[lowerIndex:upperIndex].mean()
    velComp[j] = v[lowerIndex:upperIndex].mean()
    XComp[j] = X[lowerIndex:upperIndex].mean()
    TComp[j] = T[lowerIndex:upperIndex].mean()
    csComp[j] = cs[lowerIndex:upperIndex].mean()
    '''
    rhoComp[j] = 0.5*(rho[compIndex-1]+rho[compIndex])
    pressComp[j] = 0.5*(press[compIndex-1]+press[compIndex])
    velComp[j] = 0.5*(v[compIndex-1]+v[compIndex])
    XComp[j] = 0.5*(X[compIndex-1]+X[compIndex])
    TComp[j] = 0.5*(T[compIndex-1]+T[compIndex])
    csComp[j] = 0.5*(cs[compIndex-1]+cs[compIndex])
    '''

    '''   
    plt.clf()
    fig = plt.figure(figsize=(14,9))

    plt.subplot(2,3,1)
    plt.scatter(r/pc, rho/rhoISM, s=1)
    #plt.axis([-0.01-boxSize,0.01+boxSize,0.0,rhoMax])
    # plt.axvline( x=mirrorLeft, c='k' )
    # plt.axvline( x=mirrorRight, c='k' )
    plt.xlim(0.0,xmax)
    plt.axvline(rComp, c='k')
    plt.xlabel(r'$r$ (pc)', fontsize=15)
    plt.ylabel(r'$\rho/\rho_{\rm ISM}$', fontsize=15)
    plt.title(r'$t=$'+str(time[j]/yr)+' yr')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    plt.subplot(2,3,2)
    plt.scatter(r/pc, press/pressISM, s=1)
    #plt.axis([-0.01-boxSize,0.01+boxSize,0.0,tempMax])
    # plt.axvline( x=mirrorLeft, c='k' )
    # plt.axvline( x=mirrorRight, c='k' )
    plt.xlim(0.0,xmax)
    plt.axvline(rComp, c='k')
    plt.xlabel(r'$r$ (pc)', fontsize=15)
    plt.ylabel(r'$P/P_{\rm ISM}$', fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    plt.subplot(2,3,3)
    plt.scatter(r/pc, v/cmtokm, s=1)
    #plt.axis([-0.01-boxSize,0.01+boxSize,vxMin,vxMax])
    # plt.axvline( x=mirrorLeft, c='k' )
    # plt.axvline( x=mirrorRight, c='k' )
    plt.xlim(0.0,xmax)
    plt.axvline(rComp, c='k')
    plt.xlabel(r'$r$ (pc)', fontsize=15)
    plt.ylabel(r'$v_{r}$ (km/s)', fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    plt.subplot(2,3,4)
    plt.scatter(r/pc, X, s=1)
    #plt.axis([-0.01-boxSize,0.01+boxSize,vxMin,vxMax])
    # plt.axvline( x=mirrorLeft, c='k' )
    # plt.axvline( x=mirrorRight, c='k' )
    plt.xlim(0.0,xmax)
    plt.axvline(rComp, c='k')
    plt.xlabel(r'$r$ (pc)', fontsize=15)
    plt.ylabel('Ejecta Fraction', fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    plt.subplot(2,3,5)
    plt.scatter(r/pc, T, s=1)
    #plt.axis([-0.01-boxSize,0.01+boxSize,vxMin,vxMax])
    # plt.axvline( x=mirrorLeft, c='k' )
    # plt.axvline( x=mirrorRight, c='k' )
    plt.xlim(0.0,xmax)
    plt.axvline(rComp, c='k')
    plt.xlabel(r'$r$ (pc)', fontsize=15)
    plt.ylabel(r'$T$ (K)', fontsize=15)
    plt.yscale('log')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    plt.subplot(2,3,6)
    plt.scatter(r/pc, cs, s=1)
    #plt.axis([-0.01-boxSize,0.01+boxSize,vxMin,vxMax])
    # plt.axvline( x=mirrorLeft, c='k' )
    # plt.axvline( x=mirrorRight, c='k' )
    plt.xlim(0.0,xmax)
    plt.axvline(rComp, c='k')
    plt.xlabel(r'$r$ (pc)', fontsize=15)
    plt.ylabel(r'$c_{s}$ (cm/s)', fontsize=15)
    plt.yscale('log')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    plt.tight_layout()

    saveas = saveasprefix + name
    fig.savefig( saveas )
    print 'saved figure',saveas
        
    plt.clf()
    fig = plt.figure(figsize=(14,9))

    plt.subplot(2,3,1)
    plt.scatter(mcoord/Msun, rho/rhoISM, s=1)
    plt.xlim(xminmc, xmaxmc)
    #plt.axis([-0.01-boxSize,0.01+boxSize,0.0,rhoMax])
    # plt.axvline( x=mirrorLeft, c='k' )
    # plt.axvline( x=mirrorRight, c='k' )
    plt.axvline(mcoordComp, c='k')
    plt.xlabel(r'$m/M_{\odot}$', fontsize=15)
    plt.ylabel(r'$\rho/\rho_{\rm ISM}$', fontsize=15)
    plt.xscale('log')
    plt.title(r'$t=$'+str(time[j]/yr)+' yr')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    plt.subplot(2,3,2)
    plt.scatter(mcoord/Msun, press/pressISM, s=1)
    plt.xlim(xminmc, xmaxmc)
    #plt.axis([-0.01-boxSize,0.01+boxSize,0.0,tempMax])
    # plt.axvline( x=mirrorLeft, c='k' )
    # plt.axvline( x=mirrorRight, c='k' )
    plt.axvline(mcoordComp, c='k')
    plt.xlabel(r'$m/M_{\odot}$', fontsize=15)
    plt.ylabel(r'$P/P_{\rm ISM}$', fontsize=15)
    plt.xscale('log')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    plt.subplot(2,3,3)
    plt.scatter(mcoord/Msun, v/cmtokm, s=1)
    plt.xlim(xminmc, xmaxmc)
    #plt.axis([-0.01-boxSize,0.01+boxSize,vxMin,vxMax])
    # plt.axvline( x=mirrorLeft, c='k' )
    # plt.axvline( x=mirrorRight, c='k' )
    plt.axvline(mcoordComp, c='k')
    plt.xlabel(r'$m/M_{\odot}$', fontsize=15)
    plt.ylabel(r'$v_{r}$ (km/s)', fontsize=15)
    plt.xscale('log')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    plt.subplot(2,3,4)
    plt.scatter(mcoord/Msun, X, s=1)
    plt.xlim(xminmc, xmaxmc)
    #plt.axis([-0.01-boxSize,0.01+boxSize,vxMin,vxMax])
    # plt.axvline( x=mirrorLeft, c='k' )
    # plt.axvline( x=mirrorRight, c='k' )
    plt.axvline(mcoordComp, c='k')
    plt.xlabel(r'$m/M_{\odot}$', fontsize=15)
    plt.ylabel('Ejecta Fraction', fontsize=15)
    plt.xscale('log')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    plt.subplot(2,3,5)
    plt.scatter(mcoord/Msun, T, s=1)
    plt.xlim(xminmc, xmaxmc)
    #plt.axis([-0.01-boxSize,0.01+boxSize,vxMin,vxMax])
    # plt.axvline( x=mirrorLeft, c='k' )
    # plt.axvline( x=mirrorRight, c='k' )
    plt.axvline(mcoordComp, c='k')
    plt.xlabel(r'$m/M_{\odot}$', fontsize=15)
    plt.ylabel(r'$T$ (K)', fontsize=15)
    plt.xscale('log')
    plt.yscale('log')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    plt.subplot(2,3,6)
    plt.scatter(mcoord/Msun, cs, s=1)
    plt.xlim(xminmc, xmaxmc)
    # plt.axvline( x=mirrorLeft, c='k' )
    # plt.axvline( x=mirrorRight, c='k' )
    plt.axvline(mcoordComp, c='k')
    plt.xlabel(r'$m/M_{\odot}$', fontsize=15)
    plt.ylabel(r'$c_{s}$ (cm/s)', fontsize=15)
    plt.xscale('log')
    plt.yscale('log')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    plt.tight_layout()

    saveas = saveasprefix + 'mc' + name
    fig.savefig( saveas )
    print 'saved figure',saveas
    '''
    j = j+1

vrel = np.absolute(velComp-vComp)
machComp = vrel/csComp
supersonic = machComp > 1.0
RAComp   = 2.0*G*mComp/vrel/vrel
etaComp  = 0.5*machComp*machComp/(machComp*machComp-1.0)*RAComp/radComp
lamb = 0.250 # lambda_c for gamma=5/3 from Bondi (1952)
RA_MT = 2.0*G*mComp/(csComp*csComp+vrel*vrel)
Mdot = np.pi*RA_MT*RA_MT*rhoComp*np.sqrt(lamb*lamb*csComp*csComp+vrel*vrel) # from Moeckel & Throop (2009)

Mint = np.zeros(len(Mdot))
Mcumulative = 0.0
for l in range(0,len(Mdot)) :
    Mcumulative = Mcumulative + Mdot[l]*dt
    Mint[l] = Mcumulative

fig = plt.figure(figsize=(13,10))

plt.subplot(3,4,1)
plt.plot(time/yr, rhoComp)
plt.xlabel(r'$t$ (yr)', fontsize=15)
plt.ylabel(r'$\rho$ (g/cm$^{3}$)', fontsize=15)
plt.yscale('log')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim(0.0,T_End/yr+t0)

plt.subplot(3,4,2)
plt.plot(time/yr, pressComp)
plt.xlabel(r'$t$ (yr)', fontsize=15)
plt.ylabel(r'$P$ (baryes)', fontsize=15)
plt.yscale('log')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim(0.0,T_End/yr+t0)

plt.subplot(3,4,3)
plt.plot(time/yr, vrel)
plt.xlabel(r'$t$ (yr)', fontsize=15)
plt.ylabel(r'$v$ (cm/s)', fontsize=15)
plt.yscale('log')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim(0.0,T_End/yr+t0)

plt.subplot(3,4,4)
plt.plot(time/yr, machComp)
plt.xlabel(r'$t$ (yr)', fontsize=15)
plt.ylabel(r'$\mathcal{M}$', fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim(0.0,T_End/yr+t0)

plt.subplot(3,4,5)
plt.plot(time/yr, RAComp/radComp)
plt.xlabel(r'$t$ (yr)', fontsize=15)
plt.ylabel(r'$R_{A}/R$', fontsize=15)
plt.yscale('log')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim(0.0,T_End/yr+t0)

plt.subplot(3,4,6)
plt.plot(time/yr, XComp)
plt.xlabel(r'$t$ (yr)', fontsize=15)
plt.ylabel('Ejecta Fraction', fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim(0.0,T_End/yr+t0)

plt.subplot(3,4,7)
plt.plot(time/yr, TComp)
plt.xlabel(r'$t$ (yr)', fontsize=15)
plt.ylabel(r'$T$ (K)', fontsize=15)
plt.yscale('log')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim(0.0,T_End/yr+t0)

plt.subplot(3,4,8)
plt.plot(time/yr, csComp)
plt.xlabel(r'$t$ (yr)', fontsize=15)
plt.ylabel(r'$c_{s}$ (cm/s)', fontsize=15)
plt.yscale('log')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim(0.0,T_End/yr+t0)

plt.subplot(3,4,9)
plt.scatter(time[supersonic]/yr, etaComp[supersonic])
plt.xlabel(r'$t$ (yr)', fontsize=15)
plt.ylabel(r'$\eta$', fontsize=15)
plt.yscale('log')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim(0.0,T_End/yr+t0)

plt.subplot(3,4,10)
plt.plot(time/yr, Mdot/Msun*yr)
plt.xlabel(r'$t$ (yr)', fontsize=15)
plt.ylabel(r'$\dot{M}$ ($M_{\odot}$/yr)', fontsize=15)
plt.yscale('log')
plt.ylim(1.0e-23,1.0e-14)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim(0.0,T_End/yr+t0)

plt.subplot(3,4,11)
plt.plot(time/yr, Mint/Msun)
plt.xlabel(r'$t$ (yr)', fontsize=15)
plt.ylabel(r'Accreted Mass ($M_{\odot}$)', fontsize=15)
plt.yscale('log')
plt.ylim(1.0e-13,1.0e-12)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim(0.0,T_End/yr+t0)

plt.tight_layout()

saveas = 'comptraj.pdf'
fig.savefig( saveas )
print 'saved figure',saveas

