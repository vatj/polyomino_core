import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from itertools import product

def VisualiseSingleShape(shape,ax=None,corner=(0,1),extent=1):
    cols=['darkgreen','royalblue','firebrick','goldenrod','mediumorchid']
    hatchs=['//','.','\\\\','O']
    ar_offsets={0:(0,-.25,0,.25),1:(-.25,0,.25,0),2:(0,.25,0,-.25),3:(.25,0,-.25,0)}
    rescale=False
    if ax is None:
        rescale=True
        plt.figure()
        ax = plt.gca()

    dx=shape[0]
    dy=shape[1]
    max_d=max(dx,dy)
    extentS=extent/max_d
    for i,j in product(range(dx),range(dy)):
        if(shape[2+i+j*dx]):
            
            coords=corner[0]+(i/max_d)*extent, corner[1]-((j+1)/max_d)*extent

            ax.add_patch(Rectangle(coords, extentS, extentS, facecolor=cols[(shape[2+i+j*dx]-1)//4],edgecolor='darkgray',fill=True,hatch=hatchs[(shape[2+i+j*dx]-1)//4],lw=0,transform=ax.transAxes))
            ax.add_patch(Rectangle(coords, extentS,extentS ,edgecolor='gainsboro',fill=False,lw=1.5,transform=ax.transAxes))
            theta=(shape[2+i+j*dx]-1)%4;
            #ax.arrow(coords[0]+extent/2+ar_offsets[theta][0], coords[1]+extentS/2+ar_offsets[theta][1]*extent, ar_offsets[theta][2]*extent, ar_offsets[theta][3]*extent, head_width=0.05*extent, head_length=0.1*extent, fc='k', ec='k',transform=ax.transAxes)

    if rescale:
        maxB=max(dx,dy)
        #ax.set_xlim([-maxB/2-0.5,maxB/2+0.5])
        #ax.set_ylim([-maxB/2-0.5,maxB/2+0.5])
        ax.set_aspect('equal')
        ax.set_axis_off()
        plt.show(block=False)
