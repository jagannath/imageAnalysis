






######### LINE PLOTTING FUNCTION ######################

def plot_LineIntensity(circles,imageDir,subDir):
    # Note = For accesing the images, the coordinates are read as (y,x) [rows
    # and cols]
    def __plotCircle(cChoice,cCoords):
        cCoords = circles[0][cChoice]
        print cChoice,nbrCircles
        profile = PlotLineProfile(imageDir,subDir,nbrCircles)
        dyeImg, cdyeImg = profile._openImages(profile._getImageName(4),-1)
        #axLeft,axRight = profile.axTup[cChoice],profile.axTup[cChoice+1]
        x0,y0,r =  cCoords
        print "Circle at (%f,%f); radius = %f"%(x0,y0,r)
        
        xVal,yVal,xyCoords,start,end =  _getIntensityList(dyeImg,'Xaxis',(x0,y0,r))
#        profile._overLayLine(cdyeImg,start,end,col=colors[0][0])        
#        axLeft = axLeft.imshow(cdyeImg,cmap='gray')
#        profile._plotLine(axRight,xVal,yVal,col=colors[0][1])

        #for i,ax in enumerate(profile.axTup[:-1]):
         #   xVal,yVal,xyCoords,start,end = profile.getIntensityList(dyeImg,orientation[i],(x0,y0,r))
          #  profile._overLayLine(cdyeImg,start,end,col=colors[i][0])
           # profile._plotLine(ax,xVal,yVal,col=colors[i][1])
           # profile._plotLabels(ax,orientation[i] + ' Profile','Range (Pixels)', 'Raw Intensity')
        return (cdyeImg,(xVal,yVal))
    # End of function

    def _plotLine(ax,x,y,col='#88CCEE'):
        ax.plot(x,y,linewidth=2.0,color=col)
        return ax

    def _plotLabels(ax,title,xlabel,ylabel):
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.text( 0.5,0.9,title,
                horizontalalignment='center', transform = ax.transAxes,
                size='large',weight='bold')
        return ax

    def _showFig():
        plt.tight_layout(pad=0.8,w_pad=0.8,h_pad=1.2)
        plt.show()

    def _getIntensityList(img,orientation,(x0,y0,r)):
        extfactor = 1.2
        extRad = r*extfactor
        def __checkEdges(x):
            if x > 512: x = 512
            elif x < 0: x = 0
            else: x
            return x

        if orientation is 'Xaxis':
            x1 = x0 - extRad 
            x2 = x0 + extRad
            x1,x2 = map(__checkEdges,[x1,x2])
            yx = [(int(y0),int(i)) for i in np.arange(x1,x2)]
            xaxis = np.arange(x1,x2)
            start = (int(x1),int(y0))
            end = (int(x2),int(y0))
        elif orientation is 'Yaxis':
            y1 = y0 - extRad
            y2 = y0 + extRad
            y1,y2 = map(__checkEdges,[y1,y2])
            yx = [(j,x0) for j in np.arange(y1,y2)]
            xaxis = np.arange(y1,y2)
            start = (int(x0),int(y1))
            end = (int(x0),int(y2))
        elif orientation is 'Slant':
            x1 = x0 - extRad*np.sin(np.pi/4)
            y1 = y0 - extRad*np.sin(np.pi/4)
            x2 = x0 + extRad*np.sin(np.pi/4)
            y2 = y0 + extRad*np.sin(np.pi/4)
            xN = [__checkEdges(x) for x in np.arange(x1,x2)]
            yN = [__checkEdges(y) for y in np.arange(y1,y2)]
            yx = zip(yN,xN)
            start = (int(x1),int(y1))
            end = (int(x2),int(y2))
            xaxis = [item[1] for item in yx]
        intList = [img[i] for i in yx]
        
        return (xaxis,intList,yx,start,end)


    ####### END OF INTERNAL FUNCTION #####


    colBlue = [( 51,34,136),'#332288']
    colGreen = [( 153,153,51),'#999933']
    colRed = [( 136,34,85),'#882255']
    colors = [colBlue,colGreen,colRed]
    orientation = ['Slant','Xaxis','Yaxis']
    
    nbrCircles = len(circles[0])
    fig, axTupTup =  plt.subplots(nrows=nbrCircles,ncols=2) 
    axTup = [ax for axTup in axTupTup for ax in axTup]

    for i,cCoords in enumerate(circles[0]):
        axL,axR = [],[]
        print "Processing %d ..."%(i)
        cDyeImg,(xVal,yVal) = __plotCircle(i,cCoords)
        axL,axR = axTup[i*2],axTup[(i*2)+1]
        
        axR = _plotLine(axR,xVal,yVal)
        axL = axL.imshow(cDyeImg,cmap='gray')
        axR = _plotLine(axR,xVal,yVal)
        axR = _plotLabels(axR,'Circle #'+str(i),'Raw Pixels','Raw Intensity')

    fname = os.path.join(imageDir,subDir+'_LINE_RAW_.png')
    print fname
    plt.savefig(fname)
    sys.exit(1)
    #profile._showFig()


######## RADIAL PLOTTING FUNCTION ###################
