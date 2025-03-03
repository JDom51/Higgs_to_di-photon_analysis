import RunAnalysis
import ROOT as r
from ROOT import *
import os
import sys
import math
import numpy as np
sys.path.insert(0, "backend") # allow code to be imported from subdirectory
from ShowHistogram import plotHistogram,histogramsFromFile
from operator import itemgetter


def landau(x, c):
    return c[0]*r.TMath.Landau(x, c[1], c[2])

def gaussian(x, c):
    return c[0]*np.exp(-0.5*((x-c[1])/c[2])**2)

def exp(x, c):
    return math.e**(c[0] + c[1]*x)

def poln(x, c):
    result = 0
    for n, coef in enumerate(c):
        result += x**n *coef
    return result

def load_parameters(path, name):
    params = []
    parerrors = []
    with open(path + name + ".txt", "r") as file:
        string = file.read()
    string_cheese = string.split(",\n")
    for i in string_cheese[0].split(","):
        params.append(float(i))
    for i in string_cheese[1].split(","):
        parerrors.append(float(i))
    chi2 = float(string_cheese[2].split(",")[0])
    Ndf = float(string_cheese[2].split(",")[1])
    chi2r = float(string_cheese[2].split(",")[2])
    return params, parerrors, chi2, Ndf, chi2r

def saving_fit(key, additional_identifier, fitting_params, fitting_parerrors, chi2r, fit_data):
    try:
        new_directory = f"Parameters/{additional_identifier}"
        os.mkdir(new_directory)
    except FileExistsError:
        pass
    with open(new_directory+ "/" + key + additional_identifier + ".txt", "w") as file:
        string = ""
        for i in fitting_params:
            string += str(i) + ","
        string += "\n"
        for i in fitting_parerrors:
            string += str(i) + ","
        string += "\n"
        
        string += str(fit_data.Chi2()) + "," + str(fit_data.Ndf()) + "," + str(chi2r)
        file.write(string)
        

def getting_fit_param(key, additional_identifier, minx, maxx, CURRENT_CUT):
    hist_file_atlas = r.TFile.Open("out/yy.root","READ")
    histograms_atlas, printString_atlas = histogramsFromFile(hist_file_atlas)
    fitting_params = []
    fitting_parerrors = []
    canvas = {}

    fitting_fun = {"expo": 2, "landau": 3, "pol2": 3, "pol3": 4, "pol4": 5, "pol5": 6, "pol6": 7, "pol7": 8, "pol8": 9, "pol9": 10, "gaus": 3}

    fit_data = histograms_atlas[CURRENT_CUT].Fit(key, "S", "", minx, maxx)
    histograms_atlas[CURRENT_CUT]
    canvas0 = plotHistogram(histograms_atlas[CURRENT_CUT])
    
    save_directory = f"InvMassFits/{additional_identifier}"
    try:
        os.mkdir(save_directory)
    except FileExistsError:
        pass
    
    canvas0.Print(save_directory + "/InvMassFit_" + key + additional_identifier + ".pdf")
    hist_file_atlas.Close()
    try:
        chi2r = fit_data.Chi2()/fit_data.Ndf()
    except ZeroDivisionError:
        chi2r = -1
    for i in range(0, fitting_fun[key]):
        fitting_params.append(fit_data.Parameter(i))
        fitting_parerrors.append(fit_data.ParError(i))        
    
    saving_fit(key, additional_identifier, fitting_params, fitting_parerrors, chi2r, fit_data)
    
    return fitting_params, fitting_parerrors, chi2r, canvas0


def residuals(fitting_fun, additional_identifier, fit_params, fit_paramerror, minx, maxx, CURRENT_CUT):    
    
    fitting_funtofun = {"expo": exp, "gaus": gaussian, "landau": landau, "pol2": poln, "pol3": poln, "pol4": poln, "pol5": poln, "pol6":poln, "pol7": poln, "pol8": poln, "pol9": poln}
    
    hist_file_atlas = r.TFile.Open("out/yy.root","READ")
    histograms_atlas, printString_atlas = histogramsFromFile(hist_file_atlas)

    at_bin_number = histograms_atlas[CURRENT_CUT].GetNbinsX()
    tot_bin_num = math.floor(at_bin_number*(1 - minx/maxx))
    res = TH1F("residuals", "residuals", tot_bin_num , minx, maxx)

    for i in range(1, tot_bin_num+1):
        # for i in range 0, n (where the max i value is n-1) will set bins from 1 to n (inclusive)
        w_hist = histograms_atlas[CURRENT_CUT]
        first_bin = round(at_bin_number * minx/maxx)
        at_v = w_hist.GetBinContent(first_bin + i)
        back_v = fitting_funtofun[fitting_fun](w_hist.GetBinCenter(first_bin + i), fit_params)
        res.SetBinContent(i, at_v - back_v)
    
    canvas1 = TCanvas()
    res.Draw("*H")
    save_directory = f"Residuals/{additional_identifier}"
    try:
        os.mkdir(save_directory)
    except FileExistsError:
        pass
    canvas1.Print(save_directory + "/residuals_" + fitting_fun + additional_identifier + ".pdf")
    
    hist_file_atlas.Close()
    return canvas1

def main_param(fit_fun, CURRENT_CUT, xmin, restricted):
    hist_file_atlas = r.TFile.Open("out/yy.root","READ")
    histograms_atlas, printString_atlas = histogramsFromFile(hist_file_atlas)
    h = histograms_atlas[CURRENT_CUT]
    bin_number = h.GetNbinsX()
    xmax = h.GetBinCenter(bin_number) + h.GetBinWidth(bin_number)/2
    
    mode = ["g", "l"] # 0 = getting fit params, 1 = load fit param 
    #user options
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    yyinvmin, yyinvmax = 120e3, 130e3
    selection = 1
    if restricted:
        selection = 0
    xmin = xmin
    cuts = CURRENT_CUT.split("restricted")[0]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    additional_identifier = f"_{xmin}-{xmax}_{yyinvmin}-{yyinvmax}_{bin_number}_{CURRENT_CUT}{mode[0]}" #a-b_c-d_e: a-b is our fit range, c-d is the range we have hidden, e is which cuts we are using (1st or 2nd)
    canvas0 = 0
    if not selection:
        fp, fpe, chi2r, canvas0 = getting_fit_param(fit_fun, additional_identifier, xmin, xmax, CURRENT_CUT)
    else:
        fp, fpe, chi2, Ndf, chi2r = load_parameters("Parameters/"+additional_identifier+ "local"+"/", fit_fun+additional_identifier )
    canvas1 = residuals(fit_fun, additional_identifier + mode[1], fp, fpe, xmin, xmax, CURRENT_CUT)


def get_chi2r_from_file(path, name):
    with open(path + name, "r") as file:
        string = file.read()
    string_cheese = string.split(",\n")
    return float(string_cheese[2].split(",")[2])

def ideal_chi2r(path, name):
    x = np.linspace(0,200,201)
    max_chi2 = -1
    with open(path + name, "r") as file:
        string = file.read()
    string_cheese = string.split(",\n")
    ndf = float(string_cheese[2].split(",")[1])
    maxpchi2 = 0
    for i in x:
        pchi2 = r.Math.chisquared_pdf(i, ndf)
        if pchi2 > maxpchi2:
            max_chi2 = i
            maxpchi2= pchi2
    try:
        max_chi2r = max_chi2/ndf
    except:
        max_chi2r = -1
    return max_chi2r

def get_chi2cdf_upperbound(path, name, cutoff):
    x = np.linspace(0,200,401)
    
    with open(path + name, "r") as file:
        string = file.read()
    string_cheese = string.split(",\n")
    ndf = float(string_cheese[2].split(",")[1])
    for i in x:
        cpchi2 = r.Math.chisquared_cdf_c(i, ndf)
        if cpchi2 <= cutoff:# this is as we are starting at p = 1 and going down
            return i/ndf
        
def get_chi2cdf_lowerbound(path, name, cutoff):
    x = np.linspace(0,200,401)
    
    with open(path + name, "r") as file:
        string = file.read()
    string_cheese = string.split(",\n")
    ndf = float(string_cheese[2].split(",")[1])
    for i in x:
        cpchi2 = r.Math.chisquared_cdf(i, ndf)
        if cpchi2 >= cutoff:# this is as we are starting at p = 1 and going down
            return i/ndf, ndf


def clean_cdf_checked_data(data):
    clean_data = []
    for row in data:
        if row[-1] == True:
            clean_data.append(row)
    return clean_data


def sort_file_data(data, index):
     return sorted(data, key=itemgetter(index))


def getting_file_data(path, file_names):
    data = []
    cutoff = 0.025
    for name in file_names:
        strung = name.split("/")[1].split("_")
        try:
            chi2r = get_chi2r_from_file(path, name)
            if chi2r != 0:
                chi2r_upper = get_chi2cdf_upperbound(path, name, cutoff)
                chi2r_lower, ndf = get_chi2cdf_lowerbound(path, name, cutoff)
                ichi2r = ideal_chi2r(path, name)
                if ichi2r == 0:
                    ichi2r = -1
                deviation = abs(round((chi2r-ichi2r)/ichi2r*100))
                
                passed = chi2r_lower < chi2r < chi2r_upper
                data.append([strung[0], strung[1], chi2r, deviation, chi2r_upper, name, passed])
        except IsADirectoryError:
            print("hello", name)
    return data

def get_save_data(path, file_name):
    with open(path + file_name, "r") as file:
        save_data = file.read()
    return save_data

def display_bckd_fits(data):
    for i in data:
        print("fit:", i[0], "range:", i[1], "chi2r:", i[2], "deviation:  \t" + str(i[3]) + "%\tpassed5%:" + str(i[-1]))

def get_clean_save_data(save_data, clean5_data):
    clean_save_data = {}
    for row in clean5_data:
        clean_save_data[row[-2]] = save_data[row[-2]]
    return clean_save_data

def save_save_data(path, clean_save_data):
    for file_name in clean_save_data:
        with open(path + file_name.split("/")[1], "w") as file:
            file.write(clean_save_data[file_name])

def main_pass_primary(CURRENT_CUT):
    path = "Parameters/"
    directory_names = os.listdir(path)
    save_data = {}
    file_names = []
    for directory in directory_names:
        working_file_names = os.listdir(path+directory+"/")
        for working_name in working_file_names:
            if working_name.split("_")[-1].split(".")[0][:-1] == CURRENT_CUT:
                file_names.append(directory+"/"+working_name)
    data = getting_file_data(path, file_names)
    for file_name in file_names:
        save_data[file_name] = get_save_data(path, file_name)
    ds = sort_file_data(data, 3) # sorted in terms of deviation from ideal chi
    clean_5data = clean_cdf_checked_data(ds)
    clean_save_data = get_clean_save_data(save_data, clean_5data)
    save_save_data("Passed_primary/", clean_save_data)
    #display_bckd_fits(clean_5data)
    print("\n")
    print(len(clean_5data), len(data))
    return 0

def saving_fit(key, additional_identifier, fitting_params, fitting_parerrors, chi2r, fit_data):
    try:
        new_directory = f"Parameters/{additional_identifier}local"
        os.mkdir(new_directory)
    except FileExistsError:
        pass
    with open(new_directory+ "/" + key + additional_identifier + ".txt", "w") as file:
        string = ""
        for i in fitting_params:
            string += str(i) + ","
        string += "\n"
        for i in fitting_parerrors:
            string += str(i) + ","
        string += "\n"
        
        string += str(fit_data.Chi2()) + "," + str(fit_data.Ndf()) + "," + str(chi2r)
        file.write(string)


def get_chi2(data, func, params, conversion):
    x_observed = data[0]
    y_observed = data[1]
    chi2local = 0
    for i in range(0, len(x_observed)):
        y_expected = func(conversion*x_observed[i], params)
        try:
            if y_observed[i] != 0:
                #print(y_observed[i], y_expected, x_observed[i])
                chi2local += ((y_observed[i]-y_expected)**2)/ abs(y_observed[i]) # as in root y_unc = sqrt(entries = y_obs)
        except ZeroDivisionError:
            chi2local += 0
    return chi2local


def local_residuals(fitting_fun, additional_identifier, fit_params, fit_paramerror, minx, maxx, CURRENT_CUT):    
    
    fitting_funtofun = {"expo": exp, "gaus": gaussian, "landau": landau, "pol2": poln, "pol3": poln, "pol4": poln, "pol5": poln, "pol6":poln, "pol7": poln, "pol8": poln, "pol9": poln}
    
    hist_file_atlas = r.TFile.Open("out/yy.root","READ")
    histograms_atlas, printString_atlas = histogramsFromFile(hist_file_atlas)

    at_bin_number = histograms_atlas[CURRENT_CUT].GetNbinsX()
    tot_bin_num = math.floor(at_bin_number*(1 - minx/maxx))
    
    res = TH1F("residuals", "residuals", tot_bin_num , minx, maxx)

    for i in range(1, tot_bin_num+1):
        w_hist = histograms_atlas[CURRENT_CUT]
        res.SetBinContent(i, w_hist.GetBinContent(round(at_bin_number * minx/maxx) + i)-fitting_funtofun[fitting_fun]((i)*((maxx-minx)/tot_bin_num)+minx, fit_params))
    
    canvas1 = TCanvas()
    res.Draw("C*")
    save_directory = f"Residuals/{additional_identifier}local"
    try:
        os.mkdir(save_directory)
    except FileExistsError:
        pass
    canvas1.Print(save_directory + "/localresiduals_" + fitting_fun + additional_identifier + ".pdf")
    
    hist_file_atlas.Close()
    return canvas1


#                                   LOCAL CHI SQUARED
def get_chi2rlocal(fit_fun, h, local_range, yyinvmin, yyinvmax, mode, directory_name):
    fitting_fun = {"expo": exp, "landau": landau, "pol2": poln, "pol3": poln, "pol4": poln, "pol5": poln, "pol6": poln, "pol7": poln, "pol8": poln, "pol9": poln, "gaus": gaussian}
    bin_number = h.GetNbinsX()
    xmax = h.GetBinCenter(bin_number) + h.GetBinWidth(bin_number)/2
    xbinwidth = h.GetBinWidth(0)/2
    xvals = []
    yvals = []
    for b in range(1, bin_number+1):
        xvals.append(h.GetBinCenter(b))
        yvals.append(h.GetBinContent(b))

    lower_cutoff = xvals.index(local_range[0] - xbinwidth) + 1
    upper_cutoff = xvals.index(local_range[1] - xbinwidth) + 1
    xvals = xvals[lower_cutoff:upper_cutoff]
    yvals = yvals[lower_cutoff:upper_cutoff]
    fp, fpe, chi2, Ndf, chi2r = load_parameters("Passed_primary/", directory_name.split(".txt")[0])
    ndf = len(yvals) - len(fp) - yvals.count(0)
    return get_chi2([xvals, yvals], fitting_fun[fit_fun], fp, 1), ndf

def get_save_data(path, file_name):
    with open(path + file_name, "r") as file:
        save_data = file.read()
    return save_data

def get_clean_save_data(save_data, clean5_data):
    clean_save_data = {}
    for row in clean5_data:
        clean_save_data[row[-2]] = save_data[row[-2]]
    return clean_save_data

def save_save_data2(path, clean_save_data):
    
    for file_name in clean_save_data:
        with open(path + file_name, "w") as file:
            file.write(clean_save_data[file_name])

def get_chi2cdf_upperbound2(ndf, cutoff):
    x = np.linspace(0,200,401)
    for i in x:
        cpchi2 = r.Math.chisquared_cdf_c(i, ndf)
        if cpchi2 <= cutoff:# this is as we are starting at p = 1 and going down
            return i/ndf


def get_chi2cdf_lowerbound2(ndf, cutoff):
    x = np.linspace(0,200,401)
    for i in x:
        cpchi2 = r.Math.chisquared_cdf(i, ndf)
        if cpchi2 >= cutoff:# this is as we are starting at p = 1 and going down
            return i/ndf

            
def main_pass_secondary(CURRENT_CUT, xmin, xmax):
    hist_file_atlas = r.TFile.Open("out/yy.root","READ")
    histograms_atlas, printString_atlas = histogramsFromFile(hist_file_atlas)
    h = histograms_atlas[CURRENT_CUT]
    cutoff = 0.025
    mode = "g" # l = load fit param 
    ##### user options
    yyinvmin, yyinvmax = 120e3, 130e3
    fit_parameter_range = [xmin, xmax]
    local_range = [114e3, 136e3]
    #####
    path = "Passed_primary/"
    directory_names = []
    for name in os.listdir(path):
        if name.split("_")[-1].split(".")[0][:-1] == CURRENT_CUT:
            directory_names.append(name)
    #fitting_fun = {"expo": exp, "landau": landau, "pol2": poln, "pol3": poln, "pol4": poln, "pol5": poln, "pol6": poln, "pol7": poln, "pol8": poln, "pol9": poln, "gaus": gaussian}
    function_chi2rlocals = []
    for directory_name in directory_names:
        try:
            fun = directory_name.split("_")[0]
            chi2local, ndf = get_chi2rlocal(fun, h, local_range, yyinvmin, yyinvmax, mode, directory_name)
            try:
                chi2rlocal = chi2local/ndf
            except ZeroDivisionError:
                chi2rlocal = -1
            chi2rlupper = get_chi2cdf_upperbound2(ndf, cutoff)
            chi2rllower = get_chi2cdf_lowerbound2(ndf, cutoff)
            passed = chi2rllower < chi2rlocal and chi2rlocal < chi2rlupper
            function_chi2rlocals.append([fun, local_range, chi2rlocal, ndf, directory_name, passed])
        except IsADirectoryError:
            print("IsADirectoryError")
    clean_chi2rlocals = clean_cdf_checked_data(function_chi2rlocals)
    #for func in clean_chi2rlocals:
    #    print(*func)
    
    save_data = {}
    
    for file_name in directory_names:
        try:
            save_data[file_name] = get_save_data(path, file_name)
        except IsADirectoryError:
            print("IsADirectoryError")
    clean_save_data = get_clean_save_data(save_data, clean_chi2rlocals)
    save_save_data2("Passed_secondary/", clean_save_data)
    return 0


def get_current_cut(maximum_value_x, set_of_cuts):
    restricted = True
    with open("toinvestigate_cuts_and_ranges/current_xmax.txt", "w") as file:
        file.write(str(maximum_value_x))
    with open("toinvestigate_cuts_and_ranges/current_session.txt", "w") as file:
        file.write(set_of_cuts)
    with open("toinvestigate_cuts_and_ranges/current_session.txt", "r") as file: 
        CURRENT_CUT = file.read()
    if CURRENT_CUT.split("rest")[-1] != "rict":
        restricted = False
    return CURRENT_CUT, restricted

def get_background_integral(h, fp, fpe, back_range, fun):
    fitting_fun = {"expo": exp, "landau": landau, "pol2": poln, "pol3": poln, "pol4": poln, "pol5": poln, "pol6": poln, "pol7": poln, "pol8": poln, "pol9": poln, "gaus": gaussian}
    background_events = 0
    sum_unc = 0
    bin_num = h.GetNbinsX()
    xbinwidth = h.GetBinWidth(bin_num)/2
    xvals = []
    yvals = []
    y_expected = []
    dy_expected = []
    for b in range(1, bin_num+1):
        xvals.append(h.GetBinCenter(b))

    lower_fit = h.FindBin(back_range[0] + 1e3) 
    upper_fit = h.FindBin(back_range[1] - 1e3)
    xvals = xvals[lower_fit-1:upper_fit] # lowerfit-1 to include the 119e3 (from 118e3 to 132e3)

#    for x in xvals:
#        w_bin = h.FindBin(x)
#        yvals.append(h.GetBinContent(w_bin))
#    print("xvals = ", xvals)
#    print("fp = ", fp)
    for ind, x in enumerate(xvals):
        y_expected.append(fitting_fun[fun](x, fp))
#        print(x, "Ycomp:",yvals[ind], y_expected[ind])
#        dy_expected.append(dfitting_fun[fun](x, fp, fpe)) 
    for i, ye in enumerate(y_expected):
        background_events += ye

    return background_events

def get_cross_section(N_atlas, N_background, efficiency, luminosity):
    cross_section = (N_atlas - N_background)/(efficiency*luminosity)
    return cross_section


def main_bckd(CURRENT_CUT):
    """
    Background estimates
    """
    back_range = [120e3, 130e3]
    add_to_list_file = True
    
    hist_file_atlas = r.TFile.Open("out/yy.root","READ")
    hist_file_mc = r.TFile.Open("out/Hyy.root", "READ")
    histograms_atlas, printString_atlas = histogramsFromFile(hist_file_atlas)
    histograms_mc, printString_mc = histogramsFromFile(hist_file_mc)
    h = histograms_atlas[CURRENT_CUT]
    h2 = histograms_atlas[CURRENT_CUT]
    bin_num = h.GetNbinsX()
    max_x = h.GetBinCenter(bin_num) + h.GetBinWidth(bin_num)/2
    path = "Passed_secondary/"
    file_names = []
    for name in os.listdir(path):

        if name.split("_")[-1].split(".")[0][:-1] == CURRENT_CUT + "restrict":
            file_names.append(name)
    background_fits = []
    for name in file_names:
        try:
            fp, fpe, chi2, ndf, chi2r = load_parameters(path, name.split(".txt")[0])
            fit_range = [float(name.split("_")[1].split("-")[0]), float(name.split("_")[1].split("-")[1])]
            background_events = get_background_integral(h, fp, fpe, back_range, name.split("_")[0])
            background_fits.append(background_events)
        except IsADirectoryError:
            print("IsADirectoryError")
    avg_background = np.mean(background_fits)
    std_background = np.std(background_fits)
    print("\nFINAL BACKGROUND STUFF")
    print("Average BckGrnd =", avg_background, "\nStd on BckGrnd =", std_background, "\nSqrt BckGrnd =", np.sqrt(avg_background))
    b1 = histograms_atlas[CURRENT_CUT].FindBin(back_range[0] + 1e3)
    b2 = histograms_atlas[CURRENT_CUT].FindBin(back_range[1] - 1e3)
    eff_numerator = histograms_mc[CURRENT_CUT].Integral(b1,b2)
    print("eff numerator =", eff_numerator)
    efficiency = eff_numerator/1149.53 # DENOMINATOR IS GIVEN
    luminosity = 10.064 # GIVEN CONSTANT
    N_atlas = h.Integral(b1, b2)
    tot_integral = 0
    for i in range(b1, b2+1):
        tot_integral += h.GetBinContent(i)
    CROSS_SECTION = get_cross_section(N_atlas, avg_background, efficiency, luminosity)
    de = 0 # NEED TO FIND THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print("σ =", CROSS_SECTION, "±", std_background/(efficiency*luminosity), "(sys)", "±", np.sqrt(avg_background)/(efficiency*luminosity), "(stat) fb")
    print("σ =", CROSS_SECTION/2.29, "±",  1/2.29*std_background/(efficiency*luminosity), "(sys)", "±", 1/2.29*np.sqrt(avg_background)/(efficiency*luminosity), "(stat)", "pb\nATLAS number:", N_atlas, "Avg Background:", avg_background, "Efficiency:", efficiency)
    cross_sections = []
    for i in background_fits:
        cross_sections.append(get_cross_section(N_atlas, i, efficiency, luminosity))
    current_cross_sections = ""
    with open("cross_sections_list.txt", "r") as file:
        current_cross_sections += file.read()
    passed = True
    for i in current_cross_sections.split("\n"):
        if i.split(":")[-1] == f"{cross_sections}":
            passed = False
    if passed:
        if CROSS_SECTION > 0:
            current_cross_sections += f"{CURRENT_CUT}:{cross_sections}\n"
        else:
                with open("bad_cuts.txt", "r") as file:
                    bad_cuts = file.read()
                with open("bad_cuts.txt", "w") as file:
                    file.write(bad_cuts + CURRENT_CUT+"\n")
        with open("cross_sections_list.txt", "w") as file: 
            file.write(current_cross_sections)
            
    
    with open("cross_sections_list.txt", "r") as file:
        string_working_cross_sections = file.read().split("\n")[:-1]
    working_cross_sections = []
    for i in string_working_cross_sections:
        for n, j in enumerate(i.split(":[")[1].split(",")):
            if n != len(i.split(":[")[1].split(","))-1:
                working_cross_sections.append(float(j))
            else:
                try:
                    working_cross_sections.append(float(j[:-1]))
                except ValueError:
                    working_cross_sections = working_cross_sections
    canvas0 = TCanvas()
    sigma_plot = TH1F("cross_sections", "cross_sections", 100, min(working_cross_sections)-1, max(working_cross_sections)+1)
    for cs in working_cross_sections:    
        sigma_plot.AddBinContent(sigma_plot.FindBin(cs))
    
    sigma_plot.Draw("")
    canvas0.Print("SIGMA.pdf")
    return  0

def plot_cross_sections_graph():
    with open("cross_sections_list.txt", "r") as file:
        string_working_cross_sections = file.read().split("\n")[:-1]
    working_cross_sections = []
    for i in string_working_cross_sections:
        for n, j in enumerate(i.split(":[")[1].split(",")):
            if n != len(i.split(":[")[1].split(","))-1:
                working_cross_sections.append(float(j))
            else:
                try:
                    working_cross_sections.append(float(j[:-1]))
                except ValueError:
                    working_cross_sections = working_cross_sections
    canvas0 = TCanvas()
    sigma_plot = TH1F("cross_sections", "cross_sections", 50, min(working_cross_sections)-1, max(working_cross_sections)+1)
    for cs in working_cross_sections:    
        sigma_plot.AddBinContent(sigma_plot.FindBin(cs))
    
    sigma_plot.Draw("")
    sigma_plot.Fit("gaus", "S")
    canvas0.Print("SIGMA.pdf")
        
def main():
    optimising = False
    fitting_fun = {"expo": 2, "landau": 3, "pol2": 3, "pol3": 4, "pol4": 5, "pol5": 6, "pol6": 7, "pol7": 8, "pol8": 9, "pol9": 10, "gaus": 3}
    with open("toinvestigate_cuts_and_ranges/viable_cuts.txt","r") as file:
        string_of_cuts = file.read().split("\n")[:-1]
    set_of_cuts = []
    xminima = [106e3, 108e3, 110e3]
    max_x_values = [160e3, 150e3]
    for i in string_of_cuts:
        set_of_cuts.append(i.split(":")[0])
    for x in max_x_values:
        for cut_index in range(0, len(set_of_cuts)):
            CURRENT_CUT, restricted = get_current_cut(x, set_of_cuts[len(set_of_cuts)-cut_index-1])
            CURRENT_CUT=CURRENT_CUT
            print(f"starting {CURRENT_CUT}")
            canvases = RunAnalysis.main()
            if CURRENT_CUT.split("rest")[-1] == "rict":
                for fun in fitting_fun:
                    for xmin in xminima:
                        main_param(fun, CURRENT_CUT, xmin, restricted)
                print("start_primary")
                main_pass_primary(CURRENT_CUT)
                print("start_secondary")
                for xmin in xminima:
                    for xmax in max_x_values:
                        main_pass_secondary(CURRENT_CUT, xmin, xmax)
        
            if CURRENT_CUT.split("rest")[-1] != "rict":
                main_bckd(CURRENT_CUT)
    return 0


main()
#plot_cross_sections_graph()
#canvases = None   #uncomment to hide autoplotting
