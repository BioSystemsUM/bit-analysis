
import os, shutil
from subprocess import Popen

from Bio import SeqIO


def run(organism, params, protein_file):
    inpath = './workerSubmissions/Inputs/' + organism + "/workerSubmissions/"
    workerpath = './workerSubmissions/'
    resultsWorkerPath = './resultsWorker/'
    fasta_file = inpath + "protein.faa"

    parse_protein_file(inpath, protein_file)
    create_idGenes(inpath)
    create_params(inpath, params)

    for item in os.listdir(inpath):
        s = os.path.join(inpath, item)
        d = os.path.join(workerpath, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, )
        else:
            shutil.copy2(s, d)
    #
    command = "java -jar -Xmx5G -Xss512M -XX:+HeapDumpOnOutOfMemoryError -Djavax.xml.accessExternalDTD=all -Dworkdir=workdir/ Blast.jar " + str(params["option"])
    # # #
    child = Popen(command.split(" "))

    exit_code = child.wait()

    if exit_code != 0:
        print(child.stdout)

    else:
    # # #
        if not os.path.exists('workerSubmissions/Outputs/' + organism):
            os.makedirs('workerSubmissions/Outputs/' + organism)

        for item in os.listdir("./"):
            if item == 'model.xml':
                s = os.path.join("./", item)
                d = os.path.join('workerSubmissions/Outputs/' + organism, item)
                if os.path.isdir(s):
                    shutil.copytree(s, d, )
                else:
                    shutil.move(s, d)

        dictionary = {"1":'all', "2":"selected", "3": "selected", "4": "random"}
        if 'random_number' in params.keys():
            dictionary["4"] += str(params["random_number"])
            dictionary["2"] = "random" + str(params["random_number"])

        if "additional_name" not in params.keys():
            params["additional_name"] = ''

        output_model_name =  'workerSubmissions/Outputs/' + organism + '/model_' + organism + "_" + dictionary[str(params["option"])] +"_" + params["additional_name"]
        output_model_name = output_model_name.strip("_")+'.xml'

        if os.path.exists('workerSubmissions/Outputs/' + organism + '/model.xml'):
            try:
                if os.path.exists(output_model_name):
                    os.remove(output_model_name)
                    os.rename('workerSubmissions/Outputs/' + organism + '/model.xml', output_model_name)
                else:
                    os.rename('workerSubmissions/Outputs/' + organism + '/model.xml', output_model_name)

            except Exception as e:
                print(e)

        for item in os.listdir(resultsWorkerPath):
            s = os.path.join(resultsWorkerPath, item)
            d = os.path.join('workerSubmissions/Outputs/' + organism, item)
            if os.path.isdir(s):
                shutil.copytree(s, d, )
            else:
                shutil.copy2(s, d)
            try:
                if os.path.exists('workerSubmissions/Outputs/' + organism + "/" + organism + "_" + dictionary[str(params["option"])] + "_" + params["additional_name"]+"_"  + item):
                    os.remove( 'workerSubmissions/Outputs/' + organism + "/" + organism + "_" + dictionary[str(params["option"])] + "_" + params["additional_name"]+"_"  + item)
                os.rename('workerSubmissions/Outputs/' + organism + "/" + item, 'workerSubmissions/Outputs/' + organism + "/" + organism + "_" + dictionary[str(params["option"])] + "_" + params["additional_name"]  +"_" +  item )
            except Exception as e:
                print(e)
    print("Job finished!")


def create_input_folder(organism):
    inpath = 'workerSubmissions/Inputs/' + organism + "/workerSubmissions/"
    if not os.path.exists(inpath):
        os.makedirs(inpath)


def parse_protein_file(inpath, fasta_file):
    lines = []
    print(os.getcwd())
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        lines.append(">" + str(seq_record.id) + "\n" + str(seq_record.seq) + "\n")
    lines[-1] = lines[-1].rstrip("\n")
    new_fasta_file = open(inpath + "protein.faa", "w")
    new_fasta_file.writelines(lines)
    new_fasta_file.close()


def create_idGenes(inpath):
    i = 1
    idGenes_file = open(inpath + "/idGenes" + ".txt", "w")
    lines = []
    inpath = inpath + 'protein.faa'
    for seq_record in SeqIO.parse(inpath, "fasta"):
        lines.append(str(seq_record.id) + " - " + str(i) + "\n")
        i += 1;
    lines[-1] = lines[-1].rstrip("\n")
    idGenes_file.writelines(lines)
    idGenes_file.close()


def create_params(inpath, params):
    params_file = open(inpath + "params.txt", "w")
    params_list = []
    order = ["option", "template_models", "eValue", "bitScore", "queryCoverage", "includeReactionsWithoutGPRBigg",
             "includeReactionsWithoutGPR"]
    if 'template_models' not in params.keys():
        params["template_models"] = ''

    for key in order:
        params_list.append(str(params[key]) + "\n")

    params_list[-1] = params_list[-1].rstrip("\n")
    params_file.writelines(params_list)
    params_file.close()


if __name__ == '__main__':
    directory = r"./"
    os.chdir(directory)



#######################################################
    #         STREPTOCOCCUS THERMOPHILUS     ALL
######################################################
    organism = 'Sthermophilus'
    protein_file = "proteinSthermophilus.faa"

    params = {"option": 1,
              "eValue": "1E-10", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": True, "includeReactionsWithoutGPR": True}

    create_input_folder(organism)
    run(organism, params, protein_file)

    organism = 'Sthermophilus'
    protein_file = "proteinSthermophilus.faa"

    params = {"option": 1, "additional_name": 'restrictive',
              "eValue": "1E-20", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": False, "includeReactionsWithoutGPR": False}

    create_input_folder(organism)
    run(organism, params, protein_file)
    #
    # ######################################################
    #         # STREPTOCOCCUS THERMOPHILUS     SELECTED
    # #####################################################
    organism = 'Sthermophilus'
    protein_file = "proteinSthermophilus.faa"

    params = {"option": 2,"template_models" : "iNF517;iLJ478;iYS854" ,
                 "eValue": "1E-10", "bitScore": "50", "queryCoverage": "0.75",
                  "includeReactionsWithoutGPRBigg": True, "includeReactionsWithoutGPR": True}
    #
    create_input_folder(organism)
    run(organism, params, protein_file)

    organism = 'Sthermophilus'
    protein_file = "proteinSthermophilus.faa"

    params = {"option": 2,"template_models" : "iNF517;iLJ478;iYS854" ,"additional_name": 'restrictive',
                 "eValue": "1E-20", "bitScore": "50", "queryCoverage": "0.75",
                  "includeReactionsWithoutGPRBigg": False, "includeReactionsWithoutGPR": False}
    #
    create_input_folder(organism)
    run(organism, params, protein_file)

    #######################################################
    #         STREPTOCOCCUS THERMOPHILUS     RANDOM
    ######################################################

    organism = 'Sthermophilus'
    protein_file = "proteinSthermophilus.faa"

    params = {"option": 2, 'random_number' : 1, "template_models" : "iAF987;iAM_Pc455;iNF517" ,"additional_name": 'permissive',
            "eValue": "1E-10", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": True, "includeReactionsWithoutGPR": True}

    create_input_folder(organism)
    run(organism, params, protein_file)

    organism = 'Sthermophilus'
    protein_file = "proteinSthermophilus.faa"

    params = {"option": 2, 'random_number' : 2, "template_models" : "iHN637;iAF692;iCN718" ,"additional_name": 'permissive',
            "eValue": "1E-10", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": True, "includeReactionsWithoutGPR": True}

    create_input_folder(organism)
    run(organism, params, protein_file)

    organism = 'Sthermophilus'
    protein_file = "proteinSthermophilus.faa"

    params = {"option": 2, 'random_number': 3, "template_models": "STM_v1_0;iSbBS512_1146;iYO844","additional_name": 'permissive',
              "eValue": "1E-10", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": True, "includeReactionsWithoutGPR": True}

    create_input_folder(organism)
    run(organism, params, protein_file)

    organism = 'Sthermophilus'
    protein_file = "proteinSthermophilus.faa"

    params = {"option": 2, 'random_number': 4, "template_models": "iNF517;iSDY_1059;iLJ478", "additional_name": 'permissive',
              "eValue": "1E-10", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": True, "includeReactionsWithoutGPR": True}

    create_input_folder(organism)
    run(organism, params, protein_file)

    organism = 'Sthermophilus'
    protein_file = "proteinSthermophilus.faa"

    params = {"option": 2, 'random_number': 5, "template_models": "iCN900;iIT341;iIS312_Epimastigote","additional_name": 'permissive',
              "eValue": "1E-10", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": True, "includeReactionsWithoutGPR": True}

    create_input_folder(organism)
    run(organism, params, protein_file)


    organism = 'Sthermophilus'
    protein_file = "proteinSthermophilus.faa"

    params = {"option": 2, 'random_number' : 1, "template_models" : "iAF987;iAM_Pc455;iNF517" ,"additional_name": 'restrictive',
            "eValue": "1E-20", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": False, "includeReactionsWithoutGPR": False}

    create_input_folder(organism)
    run(organism, params, protein_file)

    organism = 'Sthermophilus'
    protein_file = "proteinSthermophilus.faa"

    params = {"option": 2, 'random_number' : 2, "template_models" : "iHN637;iAF692;iCN718" ,"additional_name": 'restrictive',
            "eValue": "1E-20", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": False, "includeReactionsWithoutGPR": False}

    create_input_folder(organism)
    run(organism, params, protein_file)

    organism = 'Sthermophilus'
    protein_file = "proteinSthermophilus.faa"

    params = {"option": 2, 'random_number': 3, "template_models": "STM_v1_0;iSbBS512_1146;iYO844","additional_name": 'restrictive',
              "eValue": "1E-20", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": False, "includeReactionsWithoutGPR": False}

    create_input_folder(organism)
    run(organism, params, protein_file)

    organism = 'Sthermophilus'
    protein_file = "proteinSthermophilus.faa"

    params = {"option": 2, 'random_number': 4, "template_models": "iNF517;iSDY_1059;iLJ478", "additional_name": 'restrictive',
              "eValue": "1E-20", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": False, "includeReactionsWithoutGPR": False}

    create_input_folder(organism)
    run(organism, params, protein_file)

    organism = 'Sthermophilus'
    protein_file = "proteinSthermophilus.faa"

    params = {"option": 2, 'random_number': 5, "template_models": "iCN900;iIT341;iIS312_Epimastigote","additional_name": 'restrictive',
              "eValue": "1E-20", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": False, "includeReactionsWithoutGPR": False}

    create_input_folder(organism)
    run(organism, params, protein_file)


    #
    #######################################################
    #         XYLELLA FASTIDIOSA     ALL
    ######################################################

    organism = 'Xfastidiosa'
    protein_file = "proteinXfastidiosa.faa"

    params = {"option": 1, "eValue": "1E-10", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": True, "includeReactionsWithoutGPR": True}

    create_input_folder(organism)
    run(organism, params, protein_file)

    organism = 'Xfastidiosa'
    protein_file = "proteinXfastidiosa.faa"

    params = {"option": 1, "additional_name": 'restrictive',
              "eValue": "1E-20", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": False, "includeReactionsWithoutGPR": False}

    create_input_folder(organism)
    run(organism, params, protein_file)
    #


    #######################################################
    #         XYLELLA FASTIDIOSA     SELECTED
    ######################################################

    organism = 'Xfastidiosa'
    protein_file = "proteinXfastidiosa.faa"

    params = {"option": 2,"template_models" : "iIT341;iCN718;iJB785",
             "eValue": "1E-10", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": True, "includeReactionsWithoutGPR": True}

    create_input_folder(organism)
    run(organism, params, protein_file)

    organism = 'Xfastidiosa'
    protein_file = "proteinXfastidiosa.faa"

    params = {"option": 2,"template_models" : "iIT341;iCN718;iJB785","additional_name": 'restrictive',
             "eValue": "1E-20", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": False, "includeReactionsWithoutGPR": False}

    create_input_folder(organism)
    run(organism, params, protein_file)

    #######################################################
    #         XYLELLA FASTIDIOSA     RANDOM
    # ######################################################
    #
    organism = 'Xfastidiosa'
    protein_file = "proteinXfastidiosa.faa"

    params = {"option": 2, 'random_number': 1, "template_models": "iSSON_1240;iLB1027_lipid;iYO844",
              "additional_name": 'permissive',
              "eValue": "1E-10", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": True, "includeReactionsWithoutGPR": True}

    create_input_folder(organism)
    run(organism, params, protein_file)

    organism = 'Xfastidiosa'
    protein_file = "proteinXfastidiosa.faa"

    params = {"option": 2, 'random_number': 2, "template_models": "iYL1228;iYO844;iNF517",
              "additional_name": 'permissive',
              "eValue": "1E-10", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": True, "includeReactionsWithoutGPR": True}

    create_input_folder(organism)
    run(organism, params, protein_file)

    organism = 'Xfastidiosa'
    protein_file = "proteinXfastidiosa.faa"

    params = {"option": 2, 'random_number': 3, "template_models": "iEK1008;iAM_Pf480;iAF987",
              "additional_name": 'permissive',
              "eValue": "1E-10", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": True, "includeReactionsWithoutGPR": True}

    create_input_folder(organism)
    run(organism, params, protein_file)
    #
    organism = 'Xfastidiosa'
    protein_file = "proteinXfastidiosa.faa"

    params = {"option": 2, 'random_number': 4, "template_models": "iAM_Pf480;iYS1720;iAM_Pv461",
              "additional_name": 'permissive',
              "eValue": "1E-10", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": True, "includeReactionsWithoutGPR": True}

    create_input_folder(organism)
    run(organism, params, protein_file)
    organism = 'Xfastidiosa'
    protein_file = "proteinXfastidiosa.faa"
    #
    params = {"option": 2, 'random_number': 5, "template_models": "iAF987;iYS854;iAF692",
              "additional_name": 'permissive',
              "eValue": "1E-10", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": True, "includeReactionsWithoutGPR": True}

    create_input_folder(organism)
    run(organism, params, protein_file)

    organism = 'Xfastidiosa'
    protein_file = "proteinXfastidiosa.faa"

    params = {"option": 2, 'random_number' : 1, "template_models" : "iSSON_1240;iLB1027_lipid;iYO844" ,"additional_name": 'restrictive',
              "eValue": "1E-20", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": False, "includeReactionsWithoutGPR": False}

    create_input_folder(organism)
    run(organism, params, protein_file)

    organism = 'Xfastidiosa'
    protein_file = "proteinXfastidiosa.faa"

    params = {"option": 2, 'random_number': 2 , "template_models": "iYL1228;iYO844;iNF517","additional_name": 'restrictive',
              "eValue": "1E-20", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": False, "includeReactionsWithoutGPR": False}

    create_input_folder(organism)
    run(organism, params, protein_file)

    organism = 'Xfastidiosa'
    protein_file = "proteinXfastidiosa.faa"

    params = {"option": 2, 'random_number': 3, "template_models": "iEK1008;iAM_Pf480;iAF987","additional_name": 'restrictive',
              "eValue": "1E-20", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": False, "includeReactionsWithoutGPR": False}

    create_input_folder(organism)
    run(organism, params, protein_file)
    #
    organism = 'Xfastidiosa'
    protein_file = "proteinXfastidiosa.faa"

    params = {"option": 2, 'random_number': 4, "template_models": "iAM_Pf480;iYS1720;iAM_Pv461","additional_name": 'restrictive',
              "eValue": "1E-20", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": False, "includeReactionsWithoutGPR": False}

    create_input_folder(organism)
    run(organism, params, protein_file)
    organism = 'Xfastidiosa'
    protein_file = "proteinXfastidiosa.faa"
    #
    params = {"option": 2, 'random_number': 5, "template_models": "iAF987;iYS854;iAF692","additional_name": 'restrictive',
              "eValue": "1E-20", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": False, "includeReactionsWithoutGPR": False}

    create_input_folder(organism)
    run(organism, params, protein_file)

    #######################################################
    #         MYCOBACTERIUM TUBERCULOSIS     ALL
    ######################################################

    organism = 'Mtuberculosis'
    protein_file = "proteinMtuberculosis.faa"

    params = {"option": 1, "eValue": "1E-10", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": True, "includeReactionsWithoutGPR": True}

    create_input_folder(organism)
    run(organism, params, protein_file)

    organism = 'Mtuberculosis'
    protein_file = "proteinMtuberculosis.faa"

    params = {"option": 1, "additional_name": 'restrictive',
              "eValue": "1E-20", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": False, "includeReactionsWithoutGPR": False}

    create_input_folder(organism)
    run(organism, params, protein_file)

    # ######################################################
    # #         MYCOBACTERIUM TUBERCULOSIS     SELECTED
    # ######################################################
    #
    organism = 'Mtuberculosis'
    protein_file = "proteinMtuberculosis.faa"

    params = {"option": 2, "template_models": "iSynCJ816;iAF987;iJB785",
              "eValue": "1E-10", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": True, "includeReactionsWithoutGPR": True}

    create_input_folder(organism)
    run(organism, params, protein_file)

    organism = 'Mtuberculosis'
    protein_file = "proteinMtuberculosis.faa"
    #
    params = {"option": 2, "template_models": "iSynCJ816;iAF987;iJB785", "additional_name": 'restrictive',
              "eValue": "1E-20", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": False, "includeReactionsWithoutGPR": False}

    create_input_folder(organism)
    run(organism, params, protein_file)

    #
    # #######################################################
    # #         MYCOBACTERIUM TUBERCULOSIS     RANDOM
    # # ######################################################
    # #

    organism = 'Mtuberculosis'
    protein_file = "proteinMtuberculosis.faa"

    params = {"option": 2, 'random_number': 1, "template_models": "iCN900;iAM_Pv461;iYO844",
              "additional_name": 'permissive',
              "eValue": "1E-10", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": True, "includeReactionsWithoutGPR": True}

    create_input_folder(organism)
    run(organism, params, protein_file)

    organism = 'Mtuberculosis'
    protein_file = "proteinMtuberculosis.faa"

    params = {"option": 2, 'random_number': 2, "template_models": "iSFV_1184;iYS1720;iNF517",
              "additional_name": 'permissive',
              "eValue": "1E-10", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": True, "includeReactionsWithoutGPR": True}

    create_input_folder(organism)
    run(organism, params, protein_file)
    organism = 'Mtuberculosis'
    protein_file = "proteinMtuberculosis.faa"

    params = {"option": 2, 'random_number': 3, "template_models": "iCN718;iRC1080;iSDY_1059",
              "additional_name": 'permissive',
              "eValue": "1E-10", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": True, "includeReactionsWithoutGPR": True}

    create_input_folder(organism)
    run(organism, params, protein_file)
    organism = 'Mtuberculosis'
    protein_file = "proteinMtuberculosis.faa"

    params = {"option": 2, 'random_number': 4, "template_models": "iYL1228;iYO844;iND750",
              "additional_name": 'permissive',
              "eValue": "1E-10", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": True, "includeReactionsWithoutGPR": True}

    create_input_folder(organism)
    run(organism, params, protein_file)
    organism = 'Mtuberculosis'
    protein_file = "proteinMtuberculosis.faa"

    params = {"option": 2, 'random_number': 5, "template_models": "iSDY_1059;iYL1228;iIS312_Epimastigote",
              "additional_name": 'permissive',
              "eValue": "1E-10", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": True, "includeReactionsWithoutGPR": True}

    create_input_folder(organism)
    run(organism, params, protein_file)

    organism = 'Mtuberculosis'
    protein_file = "proteinMtuberculosis.faa"

    params = {"option": 2, 'random_number': 1, "template_models": "iCN900;iAM_Pv461;iYO844",
              "additional_name": 'restrictive',
              "eValue": "1E-20", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": False, "includeReactionsWithoutGPR": False}


    create_input_folder(organism)
    run(organism, params, protein_file)

    organism = 'Mtuberculosis'
    protein_file = "proteinMtuberculosis.faa"

    params = {"option": 2, 'random_number': 2, "template_models": "iSFV_1184;iYS1720;iNF517",
              "additional_name": 'restrictive',
              "eValue": "1E-20", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": False, "includeReactionsWithoutGPR": False}

    create_input_folder(organism)
    run(organism, params, protein_file)
    organism = 'Mtuberculosis'
    protein_file = "proteinMtuberculosis.faa"

    params = {"option": 2, 'random_number': 3, "template_models": "iCN718;iRC1080;iSDY_1059",
              "additional_name": 'restrictive',
              "eValue": "1E-20", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": False, "includeReactionsWithoutGPR": False}

    create_input_folder(organism)
    run(organism, params, protein_file)
    organism = 'Mtuberculosis'
    protein_file = "proteinMtuberculosis.faa"

    params = {"option": 2, 'random_number': 4, "template_models": "iYL1228;iYO844;iND750",
              "additional_name": 'restrictive',
              "eValue": "1E-20", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": False, "includeReactionsWithoutGPR": False}

    create_input_folder(organism)
    run(organism, params, protein_file)
    organism = 'Mtuberculosis'
    protein_file = "proteinMtuberculosis.faa"

    params = {"option": 2, 'random_number': 5, "template_models": "iSDY_1059;iYL1228;iIS312_Epimastigote",
              "additional_name": 'restrictive',
              "eValue": "1E-20", "bitScore": "50", "queryCoverage": "0.75",
              "includeReactionsWithoutGPRBigg": False, "includeReactionsWithoutGPR": False}

    create_input_folder(organism)
    run(organism, params, protein_file)