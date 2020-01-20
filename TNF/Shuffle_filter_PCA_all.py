from TNF_shuffle_sample_filter_PCA import TNF_shuffle_filter

# %%
if __name__=='__main__':
    TNF = TNF_shuffle_filter(augment_DB=True)
    TNF.read_data_dd(name='TNF_integrated_dt1e-7.h5')
    TNF.filter_data(main_RR='RR_CH4', T_min = 700, f_max=0.15, quantile=0.4, rounds = 10, PC_threshold=10)

    TNF.data_integrated_dd.compute()
    TNF.create_subset(frac=0.9)
    TNF.plot_subset(x='f_Bilger', y='RR_CH4', color_by='T')
    # TNF.plot_subset(x='f_Bilger', y='RR_H2', color_by='T')
    # TNF.plot_subset(x='f_Bilger', y='RR_CO', color_by='T')
    # TNF.plot_subset(x='f_Bilger', y='RR_OH', color_by='T')
    TNF.write_hdf(nameDB='TNF_integrated_filtered_pca',key='TNF_filtered',dt='1e-7')

    # works ... August, 2019

    plt.scatter(TNF.data_ODE_augmented_dd['f_Bilger'].compute(),TNF.data_ODE_augmented_dd['RR_CH4'].compute(),s=0.2,c=TNF.data_ODE_augmented_dd['T'].compute())
    plt.xlabel('f_Bilger')
    plt.ylabel('RR_CH4')
    plt.show()

    plt.scatter(TNF.data_ODE_augmented_dd['f_Bilger'].compute(),TNF.data_ODE_augmented_dd['T'].compute(),s=0.2)
    plt.xlabel('f_Bilger')
    plt.ylabel('T')
    plt.show()

