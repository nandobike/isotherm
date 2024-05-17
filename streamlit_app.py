# Import all necessary libraries
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from scipy.stats import linregress



st.title('Gas Adsorption Isotherm Analysis')

multi = '''Enjoy our App that integrates several gas isotherm analyses.

Questions? Suggestions? Complaints? Shoot us an email: fvb@vallejos.cl
'''
st.markdown(multi)


tab1, tab2, tab3 = st.tabs(["Load Isotherm", "BET Analysis", "SPE method"])

with tab1:
    st.header('Isotherm Data Load')
    st.markdown('Upload your isotherm as a text file. The file must only contain datapoints in ascending pressure order. Two columns separated by tabs, first for relative pressure, second for adsorbed amount in cc STP/g. See an example [here](https://raw.githubusercontent.com/nandobike/3d-vis/main/examples/a20_lao.tsv)')
    #load experimental isotherm
    #It must be a tab-separated file with two columns.
    #First column is relative pressure and second column adsorbed volume in units cc STP/g
    file = st.file_uploader("Upload isotherm file")
    if file is None:
        file = "examples/a20_lao.tsv"
        st.write(f"No file was uploaded. Loading a default isotherm file: {file}")
        #file_parameters = read_parameters(file)
        #file_name_print = file
    else:
        st.write(f'A file was uploaded: {file.name} as {file.type}')

    #if isinstance(file, str):
    exp_iso = np.genfromtxt(file, delimiter="\t") #load isotherm file into numpy array
    #print(exp_iso)

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(8,4), sharey=True, constrained_layout=True)
    #fig.tight_layout()

    for i in range(2):
        ax[i].plot(exp_iso[:,0], exp_iso[:,1],label='Experimental', marker='o', linestyle='none')
        ax[i].set_xlabel('Relative pressure P/P$_0$')
        ax[i].grid(color='aliceblue')
    ax[0].set_ylabel("Adsorbed amount (cm$^3$/g)")
    ax[0].set_ylim(bottom=0)  # adjust the bottom leaving top unchanged
    ax[1].set_xscale('log')
    ax[1].xaxis.set_major_locator(ticker.LogLocator(base=10, numticks=15))
    ax[1].set_xlim(left=1e-8, right=1.4)
    fig.suptitle('Experimental Isotherm and Interpolation to Kernel')
    st.pyplot(fig)

    #Need to add some asserts here

with tab2:
    st.header('BET Method')
    st.markdown("Let's start with the venerable BET method first. Below is the isotherm in semilog scale.")

    fig, ax = plt.subplots(figsize=(7,3))
    ax.plot(exp_iso[:,0], exp_iso[:,1],label='Experimental', marker='o', linestyle='none')
    ax.set_xlabel('Relative pressure P/P$_0$')
    ax.grid(color='aliceblue')
    ax.set_ylabel("Adsorbed amount (cm$^3$/g)")
    ax.set_ylim(bottom=0)  # adjust the bottom leaving top unchanged
    ax.set_xscale('log')
    #ax[1].xaxis.set_major_locator(ticker.LogLocator(base=10, numticks=15))
    #ax[1].set_xlim(left=1e-8, right=1.4)
    fig.suptitle('Experimental Isotherm and Interpolation to Kernel')
    st.pyplot(fig)


    st.markdown("To choose an initial range for the BET fit, let's apply the Rouquerol criterion.")
    #https://www.azom.com/article.aspx?ArticleID=10480
    fig, ax = plt.subplots(figsize=(7,3))
    ax.set_title('Rouquerol Transform')
    rouquerol = exp_iso[:,1]*(1-exp_iso[:,0])
    max_rouquerol_idx = np.argmax(rouquerol)
    ax.plot(exp_iso[:,0], rouquerol, marker='o')
    ax.plot(exp_iso[max_rouquerol_idx,0], rouquerol[max_rouquerol_idx], marker='o')
    ax.axvline(exp_iso[max_rouquerol_idx,0], color='tab:orange')
    ax.set_xscale('log')
    ax.set_xlabel("Relative pressure P/P$_0$")
    ax.grid(color='aliceblue')
    st.pyplot(fig)
   
    st.write(f'Optimum point is #{max_rouquerol_idx}, at relative pressure = {exp_iso[max_rouquerol_idx,0]}. \
             This is the recommended top limit for the BET fit.')

    sat_p = 1.00
    BET_fit_range = slice(0, max_rouquerol_idx+1)
    BET_yaxis =  exp_iso[:,0] / (exp_iso[:,1] * (sat_p - exp_iso[:,0]))
    BET_slope, BET_intercept, r, p, se = linregress(exp_iso[:,0][BET_fit_range],
                                            BET_yaxis[BET_fit_range])
    

    fig, ax = plt.subplots(figsize=(7,3))
    ax.set_title('BET plot')
    
    ax.plot(exp_iso[:,0],
            BET_yaxis, marker='o',
            label='Fit limit')
    ax.plot(exp_iso[max_rouquerol_idx,0],
            BET_yaxis[max_rouquerol_idx],
            marker='o',
            label='BET transform')
    ax.plot(exp_iso[:,0][BET_fit_range],
            exp_iso[:,0][BET_fit_range]*BET_slope+BET_intercept,
            label='BET fit')
 #   ax.axvline(exp_iso[max_rouquerol_idx,0], color='tab:orange')
 #   ax.set_xscale('log')
    ax.set_xlabel("Relative pressure P/P$_0$")
    ax.set_ylabel("BET Transform $\\frac{p}{n (p_0 - p)}$")
    ax.grid(color='aliceblue')
    ax.legend()
    ax.set_xlim(left=0, right=exp_iso[max_rouquerol_idx,0] * 1.2)
    ax.set_ylim(bottom=0, top=BET_yaxis[max_rouquerol_idx] * 1.2)
    st.pyplot(fig)

    #Probably here add an assert that intercept is positive
    #st.write(f'{BET_slope}, i={BET_intercept}')
    BET_nm = 1/(BET_slope + BET_intercept)
    st.write(f'Monolayer capacity = {BET_nm:.2f} cc STP/g')

    #STP is 0 C and 10^5 Pa https://goldbook.iupac.org/terms/view/S06036
    BET_molec_area = 0.162 #sigma, nm2
    BET_ssa = BET_nm/1e6 * 1e5 / 8.314 / 273.15 * 6.02e23 * BET_molec_area/1e9/1e9

    st.write(f'Surface area = {BET_ssa:.2f} m$^2$/g')


with tab3:

    multi = '''
    Here we use the subtracting pore effect to solve porous structure of a carbon based on our paper:

    *S. Wang, F. Vallejos-Burgos, A. Furuse, Y. Yoshikawa, H. Tanaka and K. Kaneko.* **"The subtracting pore effect method for an accurate and reliable surface area determination of porous carbons"**, Carbon, Volume 175, 2021, 77-86.

    [Here is our paper published in Carbon](https://doi.org/10.1016/j.carbon.2020.12.075)

    If you use the program please remember [to cite us](https://raw.githubusercontent.com/nandobike/isotherm/main/S0008622320312525.bib).

    Questions? Suggestions? Complaints? Shoot us an email: fvb@vallejos.cl
    '''

    st.markdown(multi)

    st.divider()
    st.header('Reference Isotherm Data')

    isotherm = './SPE_references/M-280_N2_77K.tsv'
    ref_iso = np.genfromtxt(isotherm,
                           delimiter='\t',
                           skip_header=8)

    with open(isotherm) as f:
        ref_data = []
        for i in range(6): #read the first 6 lines into a list
            ref_data.append(f.readline().strip().split('\t'))
    ref_data = dict(ref_data) #convert the list into dictionary

    for parameter in ["Temperature (K):",
                    "SSA (m2/g):",
                    "Ratio density gas/liquid:",
                    "P/P0 alpha S:"]:
        ref_data[parameter] = float(ref_data[parameter])

    col1, col2 = st.columns(2)

    with col1:
        st.write(f"Using file:  {isotherm}")

        #Show parameters in screen
        for parameter in ['Sample Name:',
                        "Gas:",
                        "Temperature (K):",
                        "SSA (m2/g):",
                        "Ratio density gas/liquid:",
                        "P/P0 alpha S:"]:
            st.text(f"{parameter} {ref_data[parameter]}")
        
        V_at_alphaS = np.interp(ref_data["P/P0 alpha S:"],
                               ref_iso[:,0],
                               ref_iso[:,1])
        st.text(f"V at alphaS: {V_at_alphaS:.2f}")

        alphaS_ref = ref_iso[:,1]/V_at_alphaS
        #print(alphaS_ref)
        

    with col2:
        fig, ax = plt.subplots(figsize=(4,3))
        ax.plot(ref_iso[:,0], ref_iso[:,1],label='Experimental', marker='o', linestyle='none')
        ax.set_xlabel('Relative pressure P/P$_0$')
        ax.grid(color='aliceblue')
        ax.set_ylabel("Adsorbed amount (cm$^3$/g)")
        ax.set_ylim(bottom=0)  # adjust the bottom leaving top unchanged
        #ax.xaxis.set_major_locator(ticker.LogLocator(base=10, numticks=15))
        #ax.set_xlim(left=1e-8, right=1.4)
        fig.suptitle('Reference Isotherm')
        st.pyplot(fig)

    st.divider()
    st.header('Alpha-S Plot')

    alphaS_sample = np.interp(exp_iso[:,0], ref_iso[:,0], alphaS_ref)
    #print(alphaS_sample)



    col1, col2 = st.columns(2)

    with col1:
        alpha_knee = st.number_input("Select Alpha value of the knee",
                                    min_value=0.0,
                                    max_value=1.5,
                                    value=0.53,
                                    step=0.01)
    with col2:
        alpha_linear = st.number_input("Select Alpha value where slope is linear",
                                    min_value=0.0,
                                    max_value=2.5,
                                    value=1.1,
                                    step=0.01)
    
    #volume at the knee
    alpha_knee_v = np.interp(alpha_knee, alphaS_sample, exp_iso[:,1])

    #fit of top curve
    alpha_fit_range = slice(np.where(alphaS_sample >= alpha_linear)[0][0],
                            alphaS_sample.shape[0])
    #print(alpha_fit_range)
    alpha_slope, alpha_intercept, r, p, se = linregress(alphaS_sample[alpha_fit_range],
                                                        exp_iso[:,1][alpha_fit_range])

    print(alpha_slope, alpha_intercept)

    fig, ax = plt.subplots(figsize=(6,4))
    ax.plot(alphaS_sample, exp_iso[:,1], label='Alpha-S', marker='o', linestyle='none')
    ax.plot([0, alpha_knee*1.8], [0, alpha_knee_v*1.8], marker=None, label='Total Area')
    ax.plot([0, alphaS_sample[-1]],
            [alpha_intercept, alphaS_sample[-1]*alpha_slope+alpha_intercept],
            label='Micropore and External')
    #ax.plot(alphaS_sample[alpha_fit_range],
    #                                                    exp_iso[:,1][alpha_fit_range])
    ax.set_xlabel('Alpha-S')
    ax.grid(color='aliceblue')
    ax.set_ylabel("Adsorbed amount (cm$^3$/g)")
    ax.set_ylim(bottom=0)  # adjust the bottom leaving top unchanged
    #ax.set_xscale('log')
    #ax.xaxis.set_major_locator(ticker.LogLocator(base=10, numticks=15))
    ax.set_xlim(left=0, right=2.5)
    ax.set_title('Alpha-S Plot')
    ax.legend()
    st.pyplot(fig)

    st.write(f'Total surface area= {alpha_knee_v/alpha_knee/V_at_alphaS*ref_data["SSA (m2/g):"]:.1f} m$^2$/g')
    st.write(f'External surface area= {alpha_slope/V_at_alphaS*ref_data["SSA (m2/g):"]:.1f} m$^2$/g')
    st.write(f'Micropore volume= {alpha_intercept/ref_data["Ratio density gas/liquid:"]:.4f} cm$^3$/g')