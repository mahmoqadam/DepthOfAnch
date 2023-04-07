import mdtraj as md
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import nglview as nv
import pandas as pd

# WORKDIR = "/Users/..."
TRAJ = "/net/orinoco/mahmoudm/betzy/aa/2e3m_replica/namd/dcd/center.dcd"
TOPOL = "/net/orinoco/mahmoudm/betzy/aa/2e3m_replica/namd/dcd/step5_charmm2namd.psf"

membraneLayerResidue = "residue 127 to 255 and not water" #2e3m
#membraneLayerResidue = "residue 129 to 256 and not water" #2e3q
proteinSelection="protein"

traj = md.load(TRAJ, top=TOPOL)

def calc_depth(traj, membraneLayerResidue, proteinSelection="protein"):

    #1. Averaging all P Z coordinate, along the trajectory
    Pplan = traj.top.select(f"({membraneLayerResidue}) and symbol P") # Get all the phosphate of the selected lipids
    average_phosphate = np.mean(traj.atom_slice(Pplan).xyz,axis=1) #calculate the average of all coordinates
    average_Z_phosphate = np.array([x[2] for x in average_phosphate]) #Average on Z only for every frames
   
    #2. Averagine the protein now.
    protSelection = traj.atom_slice(traj.top.select(f"{proteinSelection} and name CA")) #Get only the protein (AND CA)
    index = [x.residue.resSeq for x in protSelection.top.atoms] #Get the Residue number
    average_protCA = np.mean(protSelection.xyz,axis=1) #Average of all CA atoms in every frames
    average_Z_CA = np.array([x[2] for x in  average_protCA]) #Average of Z only atoms in every frames

    allz = protSelection.xyz[:,:,2]

    #nframe = protSelection.n_frames
    depth_res_time = []
    for i in range(allz.shape[1]):
        depth = average_Z_phosphate - allz[:,i]
        depth_res_time.append(depth)
    depth_res_time = np.asarray(depth_res_time)

    depth_res_time
    df = pd.DataFrame(depth_res_time, index=index)
    df = df*10

    return (df)

# to be deleted!
def calc_center_of_mass(traj):
    """Compute the center of mass by frame of trajectory compared to ref structure

    returns: returns list of tuples: (frame_index, com)
    """
    com_list = md.compute_center_of_mass(traj)
    frame_index = list(range(1, len(com_list) + 1))
    return list(zip(frame_index, com_list))

print(proteinSelection)
all_dataframe = calc_depth(traj[3000:21000], membraneLayerResidue, proteinSelection)

def average_in_time(data, title="", tick_every=10):
    dataplot = data.mean(axis=1).reset_index()
    dataplot = dataplot.rename(columns={"index":"Residue", 0:"Depth"})
    dataplot["std"] = data.std(axis=1).reset_index()[0]
    fig, ax = plt.subplots(figsize=(10,5))

    minres = dataplot.Residue.min()
    maxres = dataplot.Residue.max()
    graph = sns.lineplot(x="Residue", y="Depth", data=dataplot, color="blue", ax=ax)
    graph.set(ylabel=r"Depth of Anchoring ($\AA$)", xlabel="Residue", title=title)
    graph.fill_between(dataplot.Residue, dataplot.Depth-dataplot["std"],dataplot.Depth+dataplot["std"],color="#CCCCCC")
    hline = plt.hlines(y=0, xmin=minres, xmax = maxres, color="black", linestyles='dashed')

    graph.set_xticks(range(minres, maxres, tick_every))

    #Now the legend
    plt.legend(labels=["Average depth of anchoring","Standard deviation", "Phosphate plane"])

    dataplot.to_csv('DOA-ave.csv')



average_in_time(all_dataframe, title="Average Depth of anchoring per residue")

all_dataframe.to_csv('DOA-all.csv')

# Color background
#plt.axvspan(469, 476, alpha=0.5, color='linen')

plt.savefig("myfigure.tiff", format="tiff", transparent=False, pil_kwargs={"compression": "tiff_lzw"}, dpi=300)


plt.show()

