def Analyse(t, weighting):
    """
    this analysis is for investigating the invariant mass plot
    """
    selection_cuts = {}
    with open("toinvestigate_cuts_and_ranges/viable_cuts.txt", "r") as file:
        save_string = file.read().split("\n")[:-1]
    for name in save_string:
        selection_cuts[name.split(":")[0]] = name.split(":")[1]
    with open("toinvestigate_cuts_and_ranges/current_session.txt", "r") as file:
        name = file.read()
    with open("toinvestigate_cuts_and_ranges/current_xmax.txt", "r") as file:
        factor = float(file.read())/1000
    xmin = 0
    multiplier = 0.5
    t.SetAlias("yyInvMass","sqrt((2*photon_pt[0]*photon_pt[1])*(cosh(photon_eta[0]-photon_eta[1])-cos(photon_phi[0]-photon_phi[1])))") 
    histogram(t=t, weighting = weighting, variable = "yyInvMass", hist_id = name,   # variable thats investigated
              n_bins=multiplier*(factor-xmin), xmin=xmin *1000, xmax=factor*1000,
              cuts= f"(photon_n == 2 && photon_isTightID[1] == 1 && photon_isTightID[0] == 1 {selection_cuts[name]})")

def Analyse1(t, weighting):
    """
    for a more general analyiss function where we dont have to constantly refresh this
    """
    t.SetAlias("yyInvMass","sqrt((2*photon_pt[0]*photon_pt[1])*(cosh(photon_eta[0]-photon_eta[1])-cos(photon_phi[0]-photon_phi[1])))") 
    additional_cuts = get_additional_cuts()
    variables = get_working_var()
    a = 0
    cut_steps = 2
    #variables = {"photon_etcone20[0]": 20e3}
    with open("investigating_cuts_and_ranges/current_session.txt", "r") as file:
        current_session = int(file.read())
    saving_array = []
    for var in variables:
        additional_cuts_fixed = ""
        additional_cuts_unfixed = []
        for additional_cut in additional_cuts:
            part = additional_cut.split(":")
            if part[2] == "y":
                additional_cuts_fixed += f"&& {part[0]} {part[1]} {part[-1]}"
            else:
                additional_cuts_unfixed.append(additional_cut)
        final_w_cuts = additional_cuts_fixed
        histogram(t=t, weighting = weighting, variable = var, hist_id = var,   # variable thats investigated
              n_bins=800, xmin=a, xmax= variables[var],
              cuts = f"(photon_n == 2 && photon_isTightID[1] == 1 && photon_isTightID[0] == 1 && 120e3 <= yyInvMass && yyInvMass <= 130e3{final_w_cuts})", title = var)
        saving_array.append(var + ":" + final_w_cuts)
    saving_string = ""
    for i in saving_array:
        saving_string += i + "\n"
    with open(f"investigating_cuts_and_ranges/{currrent_session}/{current_session}.txt", "w") as file:
        file.write(saving_string)
