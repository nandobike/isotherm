# Import all necessary libraries
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker




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
        #for line in filename:
        #    st.write(line.decode('cp1252').rstrip())
        #file_parameters = read_parameters_uploaded(file)
        #file_name_print = file.name
        #1+1
        st.write(f'A file was uploaded: {file.name} as {file.type}')

    #if isinstance(file, str):
    exp_iso = np.genfromtxt(file, delimiter="\t") #load isotherm file into numpy array


    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(8,4), sharey=True, constrained_layout=True)
    #fig.tight_layout()

    for i in range(2):
        ax[i].plot(exp_iso[:,0], exp_iso[:,1],label='Experimental', marker='o', linestyle='none')
        ax[i].set_xlabel("Relative pressure P/P$_0$")
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
    st.markdown("Let's start with the venerable BET method first.")

    fig, ax = plt.subplots(figsize=(7,3))
    ax.set_title('Rouquerol Transform')
    ax.plot(exp_iso[:,0], exp_iso[:,1]*(1-exp_iso[:,0]), marker='o')
    ax.set_xscale('log')
    ax.set_xlabel("Relative pressure P/P$_0$")
    st.pyplot(fig)


with tab3:
    multi = '''This is the app that solves the microstructure of a porous carbon based on the paper:

    *S. Wang, F. Vallejos-Burgos, A. Furuse, Y. Yoshikawa, H. Tanaka and K. Kaneko et al.* **The subtracting pore effect method for an accurate and reliable surface area determination of porous carbons**, Carbon, Volume 175, 2021, 77-86.

    [Here is our paper published in Carbon](https://doi.org/10.1016/j.carbon.2020.12.075)

    If you use the program please remember [to cite us](Insert bib).

    Questions? Suggestions? Complaints? Shoot us an email: fvb@vallejos.cl
    '''
    st.markdown(multi)