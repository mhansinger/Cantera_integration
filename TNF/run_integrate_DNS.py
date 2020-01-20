from Cantera_ODE_TNF import Cantera_ODE_TNF

if __name__ == '__main__':
    myReact = Cantera_ODE_TNF()
    myReact.set_tables(name='TNF_states_DNS.h5',path='/home/hansinger/hansinger_share/TNF_database/DNS_data', key='TNF_raw_DNS')
    myReact.loop_ODE(remove_T_below=700,steps=1)
    #myReact.filter_shuffle_data(condition='RR_CH4',threshold=3)
    myReact.write_hdf(nameDB='TNF_data_integrated_DNS',key='TNF_DNS_integrated')
    # myReact.integrate_Ode(1000)

    #reactor = myReact.integrate_cantera(iloc=0)