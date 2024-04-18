#%% preparation profiles 2322 545
#I want first to print heating schemes for shots listed here, then I want to produce profiles at each step of heating scheme for each shot.
#I will use the class TCVShot to do so.


from diag_tcv.shotAnalysis.dischargeInfoMdsObject import TCVShot, plot_profiles_comparison



### Bloc ohmic heating ###
# shot=80156
# a=TCVShot(shot)
# a.plot_summary()
# plot_profiles_comparison([shot], [1.3])


### Bloc ECRH heating ###
# shot=80163
# a=TCVShot(shot)
# a.plot_summary()
# plot_profiles_comparison([shot], [1.3])

# shot=80336
# a=TCVShot(shot)
# a.plot_summary()
# plot_profiles_comparison([shot], [1])

### Bloc mixed heating ###
# shot=80162
# a=TCVShot(shot)
# a.plot_summary()
# plot_profiles_comparison([shot, shot, shot, shot], [1, 1.2, 1.6, 1.85])

# shot=80257
# a=TCVShot(shot)
# a.plot_summary()
# plot_profiles_comparison([shot, shot, shot, shot], [1, 1.2, 1.6, 1.85])

# shot=80322
# a=TCVShot(shot)
# a.plot_summary()
# plot_profiles_comparison([shot, shot, shot, shot], [1, 1.2, 1.6, 1.8])


### Bloc toroidal rotation ###

# shot=80324
# a=TCVShot(shot)
# a.plot_summary()
# plot_profiles_comparison([shot, shot], [0.9, 1.4])

# shot=80328
# a=TCVShot(shot)
# a.plot_summary()
# plot_profiles_comparison([shot, shot], [0.9, 1.4])

# shot=80745
# a=TCVShot(shot)
# a.plot_summary()
# plot_profiles_comparison([shot, shot, shot], [0.8, 1.2, 1.7])
            
# shot=80753
# a=TCVShot(shot)
# a.plot_summary()
# plot_profiles_comparison([shot, shot, shot], [0.8, 1.2, 1.7])


### Bloc density ramp ###

# shot=80376
# a=TCVShot(shot)
# a.plot_summary()
# plot_profiles_comparison([80376, 80376], [1, 1.5])