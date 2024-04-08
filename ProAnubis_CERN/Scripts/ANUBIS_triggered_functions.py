from PIL import Image
import h5py
import anubisPlotUtils as anPlot
import json
import numpy as np
import os
import hist as hi
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
matplotlib.use('TkAgg')  # or 'Qt5Agg', 'GTK3Agg', etc.
import mplhep as hep
hep.style.use([hep.style.ATLAS])
import sys
import ANUBIS_triggered_functions as ANT
import pandas as pd
import matplotlib.backends.backend_pdf
from matplotlib.ticker import MultipleLocator
import plotly
import plotly.graph_objs as go
import plotly.express as px

current_directory=  os.path.dirname(os.getcwd())

#DATA READING

def importFromTextFile(filename):
    inputText = open(filename)
    thisEvent = []
    data = [[] for tdc in range(5)]
    tdc=0
    for line in inputText:
        if "Header" in line:
            thisEvent = []
            tdc = int(line.split(" ")[1].strip(","))
        elif "Data" in line:
            thisEvent.append(int("0x"+line.split(" ")[2].strip("."),0))
        elif "EOB" in line:
            data[tdc].append(thisEvent)
    return data

def importFromHDF5File(filename):
    inputHDF5 = h5py.File(filename)
    thisEvent = []
    data = [[] for tdc in range(5)]
    tdc=0
    for event in inputHDF5['data']:
        tdc = event[0]-60928
        thisEvent = []
        for hit in event[2]:
            thisEvent.append(hit)
        data[tdc].append(thisEvent)
    return data

def importDatafile(filename):
    if "txt" in filename.split(".")[-1]:
        return importFromTextFile(filename)
    elif "h5" in filename.split(".")[-1]:
        return importFromHDF5File(filename)
    else:
        print("File type not recognized. Expect .txt or .h5 input file.")

def convertDataToHDF5(data, outfileName):
    npformattedData = []
    dt = np.dtype([('tdc', np.int_, 32), ('time', np.int_, 32),('data', np.ndarray)])
    for event in range(len(data[0])):
        for tdc in range(5):
            thisArr = np.array(data[tdc][event],dtype=np.uint32)
            thisPoint = (int(60928+tdc),int(0),thisArr)
            npformattedData.append(thisPoint)
    h5formattedData = np.array(npformattedData,dtype=dt)
    hf = h5py.File(outfileName+'.h5', 'w')
    #Read an input HDF5 file to get the right dtype. Does not work just trying to create an identical one - not sure what's missing
    preMadeDtype = h5py.File(current_directory +'\ProAnubis_CERN\ProAnubisData\240224_2158_2m_trigTest2_Triggered_RAW.h5','r')['data'].dtype
    dset = hf.create_dataset('data',(len(h5formattedData),),dtype=preMadeDtype)
    for idx in range(len(h5formattedData)):
        dset[idx] = h5formattedData[idx]
    hf.close()

def countChannels(events):
    #Expects events from one TDC, counts how many hits each channel has within the event list
    chanCounts = [0 for x in range(128)]
    for event in events:
        for word in event:
            try:
                #Extract the channel a specific hit corresponds to.
                #Every time a specific channel appears to an event, increase its hits by 1.
                chanCounts[(word>>24)&0x7f]=chanCounts[(word>>24)&0x7f]+1
            except:
                print(word>>24)
    return chanCounts

def getEventTimes(events):
    eventTimes = []
    for event in events:
        for word in event:
            #Extract (local) time of hit (UNITS?)
            eventTimes.append(word&0xfffff)
    return eventTimes

def makeSingleLayer(data, name):
    #Heatmap plot of one RPC layer. Takes already-split heat map, used by event display
    fig, ax = plt.subplots(1, figsize=(16, 8), dpi=100)
    channels= [x-0.5 for x in range(len(data)+1)]
    if(len(data)==32):
        histArr = (np.array([data]),np.array([0,1]),np.array(channels))
    else:
        histArr = ((np.array([data])).transpose(),np.array(channels),np.array([0,1]))
    thisHist = hep.hist2dplot(histArr,norm=colors.LogNorm(0.1,2))
    thisHist.cbar.remove()
    if(len(data)==32):
        plt.ylim(len(data)-0.5,-0.5)
    plt.ylabel(" ")
    plt.xlabel(" ")
    #plt.title(name)
    
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    fig.tight_layout()
    plt.savefig(current_directory+"\ProAnubis_CERN\Figures"+name+".png")
    return current_directory+"\ProAnubis_CERN\Figures"+name+".png"

def stackAndConvert(images, name="testDisplay"):
    #PIL hacking to distort images and put them together to make a primitive replication of the detector
    img = Image.open(images[0])
    total_width = 3*img.size[0]
    max_height = int(4*img.size[1])
    new_im = Image.new('RGB', (total_width, max_height))
    newData = new_im.load()
    x_offset = 0
    y_offset = 6*int(max_height/8.)
    for y in range(max_height):
        for x in range(total_width):
            #newData forms the background of the image, need to set it to all-white to start. Probably some better way to do this?
            newData[x, y] = (255, 255, 255)
    for idx, image in enumerate(images):
        img = Image.open(image)
        img = img.convert("RGBA")
        temp_im = Image.new('RGBA', (3*img.size[0], img.size[1]))
        temp_im.paste(img, (int(img.size[0]/2.),0))
        temp_im = temp_im.transform(temp_im.size, Image.AFFINE, (0.5, 1., 0, 0, 1, 0))
        pixdata = temp_im.load()
        width, height = temp_im.size
        for y in range(height):
            for x in range(width):
                if pixdata[x, y] == (255, 255, 255, 255):
                    #Manually make any white pixel transparent so that they can stack together nicely.
                    pixdata[x, y] = (255, 255, 255, 0)
        new_im.paste(temp_im, (0, y_offset), temp_im)
        y_offset = y_offset-int(max_height/28.)
        if idx == 5 or count==7:
            #Counts from the bottom up, want bigger gaps between the different chambers
            y_offset = y_offset-5*int(max_height/28.)                   
    new_im.save(current_directory+"\ProAnubis_CERN\Figures"+name+".png"+name.strip(" ")+".png", "PNG")

def makeEventDisplay(eventData,name):
    #Expects a single event, divided as [tdc0,tdc2,...,tdc4]
    countOne = countChannels([eventData[0]])
    countTwo = countChannels([eventData[1]])
    countThree = countChannels([eventData[2]])
    countFour = countChannels([eventData[3]])
    countFive = countChannels([eventData[4]])
    singEventPlots = []
    singEventPlots.append(makeSingleLayer(countOne[0:32],"Eta Triplet Low, Three Coincidences Required"))
    singEventPlots.append(makeSingleLayer(countOne[32:96],"Phi Triplet Low, Three Coincidences Required"))
    singEventPlots.append(makeSingleLayer(countOne[96:128],"Eta Triplet Mid, Three Coincidences Required"))
    singEventPlots.append(makeSingleLayer(countTwo[0:64],"Phi Triplet Mid, Three Coincidences Required"))
    singEventPlots.append(makeSingleLayer(countTwo[64:96],"Eta Triplet Top, Three Coincidences Required"))
    singEventPlots.append(makeSingleLayer(countTwo[96:128]+countThree[0:32],"Phi Triplet Top, Three Coincidences Required"))
    singEventPlots.append(makeSingleLayer(countThree[32:64],"Eta Singlet, Three Coincidences Required"))
    singEventPlots.append(makeSingleLayer(countThree[64:128],"Phi Singlet, Three Coincidences Required"))
    singEventPlots.append(makeSingleLayer(countFour[0:32],"Eta Doublet Low, Three Coincidences Required"))
    singEventPlots.append(makeSingleLayer(countFour[32:96],"Phi Doublet Low, Three Coincidences Required"))
    singEventPlots.append(makeSingleLayer(countFour[96:128],"Eta Doublet Top, Three Coincidences Required"))
    singEventPlots.append(makeSingleLayer(countFive[0:64],"Phi Doublet Top, Three Coincidences Required"))
    stackAndConvert(singEventPlots,name)
    for plot in singEventPlots:
        #Remove all the temporary plots. There's probably a better way to move pil images around than making and deleting .png files.
        os.remove(plot)

def GetEvent(eventData, num):
    return [eventData[0][num],eventData[1][num],eventData[2][num],eventData[3][num],eventData[4][num]] 

#DATA PROCESSING

def divideHitCountsByRPC(data):
    #Divides the number of hits in each channel into individual RPCs
    etaHits = [[],[],[],[],[],[]]
    phiHits = [[],[],[],[],[],[]]
    # i=0
    for event in range(0,len(data[0])):
        tdcCounts = [countChannels([data[tdc][event]]) for tdc in range(5)]
        # i+=1
        # print(i)
        etaHits[0].append(tdcCounts[0][0:32]) #Triplet Eta Low
        phiHits[0].append(tdcCounts[0][32:96]) #Triplet Phi low
        etaHits[1].append(tdcCounts[0][96:128]) #Triplet Eta Mid
        phiHits[1].append(tdcCounts[1][0:64]) #Triplet Phi Mid
        etaHits[2].append(tdcCounts[1][64:96]) #Triplet Eta Top
        phiHits[2].append(tdcCounts[1][96:128]+tdcCounts[2][0:32]) #Triplet Phi Top
        etaHits[3].append(tdcCounts[2][32:64])#Singlet Eta
        phiHits[3].append(tdcCounts[2][64:128])#Singlet Phi
        etaHits[4].append(tdcCounts[3][0:32])#Double Eta low
        phiHits[4].append(tdcCounts[3][32:96])#Double Phi Low
        etaHits[5].append(tdcCounts[3][96:128])#Doublet Eta top
        phiHits[5].append(tdcCounts[4][0:64])#Doublet Phi top
    return etaHits,phiHits

def divideEventsByRPC(data):
    #Divides the data into individual RPCs, preserving the initial words within each event
    splitEvents = []
    for event in range(len(data[0])):
        thisEvent = [[] for rpc in range(12)]
        for tdc in range(5):
            for word in data[tdc][event]:
                chan = word>>24
                if(tdc==0):
                    if(chan<32):
                        #Triplet Eta Low
                        thisEvent[0].append(word)
                    elif(chan<96):
                        #Triplet Phi Low
                        thisEvent[1].append(word)
                    else:
                        #Triplet Eta Mid
                        thisEvent[2].append(word)
                elif(tdc==1):
                    if(chan<64):
                        #Triplet Phi Mid
                        thisEvent[3].append(word)
                    elif(chan<96):
                        #Triplet Eta Top
                        thisEvent[4].append(word)
                    else:
                        #Triplet Phi 1 Top
                        thisEvent[5].append(word)
                elif(tdc==2):
                    if(chan<32):
                        #Triplet Phi 2 Top, this seems wrong too. This should be stored in index 5 not 6. Relabelled to 5
                        thisEvent[5].append(word)
                    elif(chan<64):
                        #Singlet Eta, Relabelled from 7 to 6.
                        thisEvent[6].append(word)
                    else:
                        #Singlet Phi, Relabelled from 8 to 7
                        thisEvent[7].append(word)   
                elif(tdc==3):
                    if(chan<32):
                        #Doublet Eta Low This seems wrong, why are we currently storing Doublet Eta Low and Singlet Phi in the same index 8? Relabbeld Singlet Phi index, should be fixed.
                        thisEvent[8].append(word)
                    elif(chan<96):
                        #Doublet Phi Low, This seems wrong, Doublet Phi low spans to 96 channel. Changed from 64 to 96.
                        thisEvent[9].append(word)
                    else:
                        thisEvent[10].append(word)                
                elif(tdc==4):
                    #Doublet Phi Top
                    thisEvent[11].append(word)
        splitEvents.append(thisEvent)
    return splitEvents

def countCoincidences(primRpc,SecRpc, winSize=5, minHit=3, maxHit=10):
    #Count events where at least minHit channels have an RPC hit within winSize of each channel, while the primary and secondary RPC has more than minHit and less than maxHit to filter noise events
    #Uses the pre-made hit counts instead of the event words directly
    coincArray = [0 for channel in range(len(primRpc[0]))]
    for idx, event in enumerate(primRpc):
        for channel in range(len(event)):
            nHits = 0
            for itr in range(channel-int(winSize/2.), channel+int(winSize/2.)+1):
                if itr>=0 and itr<len(event):
                    if event[channel]>0:
                        nHits = nHits+1
            if nHits>=minHit and sum(event)<=maxHit:
                if sum(SecRpc[idx])>=minHit and sum(SecRpc[idx])<=maxHit:
                    coincArray[channel]=coincArray[channel]+1
    return coincArray

def countChannels_Timed(events):
    #Expects events from one TDC, counts how many hits each channel has within the event list
    chanCounts = [[] for x in range(128)]
    for event in events:
        for word in event:

            #get time of event
            time = word&0xfffff

            try:
                #Append time of event to specific channel triggered to get event.
                chanCounts[(word>>24)&0x7f].append(time)
            except:
                print(word>>24)
    return chanCounts

def divideHitCountsByRPC_Timed(data):
    #Divides the number of hits in each channel into individual RPCs
    etaHits = [[],[],[],[],[],[]]
    phiHits = [[],[],[],[],[],[]]
    for event in range(0,len(data[0])):
        tdcCounts = [countChannels_Timed([data[tdc][event]]) for tdc in range(5)]
        etaHits[0].append(tdcCounts[0][0:32]) #Triplet Eta Low
        phiHits[0].append(tdcCounts[0][32:96]) #Triplet Phi low
        etaHits[1].append(tdcCounts[0][96:128]) #Triplet Eta Mid
        phiHits[1].append(tdcCounts[1][0:64]) #Triplet Phi Mid
        etaHits[2].append(tdcCounts[1][64:96]) #Triplet Eta Top
        phiHits[2].append(tdcCounts[1][96:128]+tdcCounts[2][0:32]) #Triplet Phi Top
        etaHits[3].append(tdcCounts[2][32:64])#Singlet Eta
        phiHits[3].append(tdcCounts[2][64:128])#Singlet Phi
        etaHits[4].append(tdcCounts[3][0:32])#Double Eta low
        phiHits[4].append(tdcCounts[3][32:96])#Double Phi Low
        etaHits[5].append(tdcCounts[3][96:128])#Doublet Eta top
        phiHits[5].append(tdcCounts[4][0:64])#Doublet Phi top
    return etaHits,phiHits

def reform_data(x,RPC,side):
    #Reformatting hit data in a specific way for coincident function to work.
    out = []
    for ind,n in enumerate(x):
        #n is list of hit times for a specific channel in an event.
        if n:
            for y in n:
                #y is hit time in hittimes n
                out.append([RPC,ind,y,side])
    return out

def coincidence(time,event_sorted,time_window):
    #Find all hits that are coincident.
    hit_locations = []
    t = time

    #Elements in event_sorted look like [RPC,CHANNEL,HIT_TIME,side]
    for hit in event_sorted:
        if t<= hit[2] < t+time_window:
            hit_locations.append(hit)
            t= hit[2] 
            #Update time for filter, the cluster will continue if more hits are coincident after this event within the time window.
    
    final_time = t+ time_window

    return hit_locations, final_time

def FindCoincidentHits(etaHits,phiHits,time_window):

    #time_window in nanoseconds

    i = 0
    coincident_hits = []

    for event in range(0,len(etaHits[0])):
        
        #Progress Bar
        i+=1
        print(f"{i/len(etaHits[0])*100:.0f}%")

        #Scan over all events in the data collection

        channels = []

        for RPC in range(6):
            #Scan through hits of each RPC

            hits_eta = reform_data(etaHits[RPC][event],RPC,'eta') #Extracting relevant element of etaHits
            hits_phi = reform_data(phiHits[RPC][event],RPC,'phi') #Extracting relevant element of phiHits
            
            #Combine all data from different RPCs in an event
            channels += hits_eta
            channels += hits_phi

        #Sort all hits in an event by time in ascending order. 
        event_sorted = sorted(channels, key=lambda x: x[2])

        #event_sorted = [[RPC,CHANNEL,HITTIME],[...],...] , with HITTIME SORTED IN ASCENDING ORDER!

        #initiliase time to time of first event.

        t = event_sorted[0][2]
        while t < 1250:
            for hit in event_sorted:
                if hit[2] == t:
                    hit_locations, final_time = coincidence(t,event_sorted,time_window)
                    coincident_hits.append([f'Event {i}',t, hit_locations]) # t here is the initial time of the cluster.
                    t = final_time# Update time 
                    
            t+= 1

    return coincident_hits

def dark_clustering(dark_coincidence, anomalous_cutoff =10):

    phi_cluster_distribution = [[] for _ in range(6)] # [[RPC1],[RPC2],[RPC3],...]
    eta_cluster_distribution = [[] for _ in range(6)]
    anomalous_clusters = [[] for _ in range(6)] #Store anomalous hits for each RPC

    all_rpc_phi_clusters = []
    all_rpc_eta_clusters = []

    for coincidence_event in dark_coincidence:

        # coincidence : ['Event x', TIMEBIN, [hit_locations]]
        hit_locations = coincidence_event[2]
        #hit_locations = [[RPC,CHANNEL,HIT_TIME,eta/phi]...]

        #Extract hit_locations in phi and eta directions.
        phi_locations = [x for x in hit_locations if x[3]=='phi']
        eta_locations = [x for x in hit_locations if x[3]=='eta']

        #Sort by channels
        phi_locations = sorted(phi_locations, key=lambda x: x[1])
        eta_locations = sorted(eta_locations, key=lambda x: x[1])


        for RPC in range(6):
            #Work out the cluster distribution for each RPC during this specific dark count event.

            rpc_phi_clusters = []
            rpc_eta_clusters = [] 

            i = 0
            for index,hit in enumerate([x for x in phi_locations if x[0]==RPC]):
                if index==0:
                    previous_element = hit[1]
                    rpc_phi_clusters.append([hit])
                else:
                    if abs(hit[1] - previous_element) > 1:
                        # Hit is not part of the same cluster, intiate a new cluster
                        rpc_phi_clusters.append([hit])
                        i += 1
                    else:
                        # Hit is part of the same cluster
                        rpc_phi_clusters[i].append(hit)
                    previous_element = hit[1]

            j = 0
            for index,hit in enumerate([x for x in eta_locations if x[0]==RPC]):
            
                if index == 0:
                    previous_element = hit[1]
                    rpc_eta_clusters.append([hit])
                else:
                    if abs(hit[1] - previous_element) > 1:
                        # Hit is not part of the same cluster
                        rpc_eta_clusters.append([hit])
                        j += 1
                    else:
                        # Hit is part of the same cluster
                        rpc_eta_clusters[j].append(hit)
                    previous_element = hit[1]

            if rpc_phi_clusters:
                all_rpc_phi_clusters.append(rpc_phi_clusters)
            if rpc_eta_clusters:
                all_rpc_eta_clusters.append(rpc_eta_clusters)

            for x in rpc_phi_clusters:
                if len(x) > anomalous_cutoff:
                    #ANOMALOUS EVENT, NOISE BURST? 
                    anomalous_clusters[RPC].append(x)
                else:
                    phi_cluster_distribution[RPC].append(len(x))
                
            for y in rpc_eta_clusters:
                if len(y) > anomalous_cutoff:
                    anomalous_clusters[RPC].append(y)
                else:
                    eta_cluster_distribution[RPC].append(len(y))

    return all_rpc_eta_clusters,all_rpc_phi_clusters,phi_cluster_distribution,eta_cluster_distribution,anomalous_clusters

def cluster(coincident_hits):
#Entries into coincident_hits have already undergone temporal clustering.
    # Now will undergo spatial clustering.

    coincident_hits_clustered = []

    for coincidence_event in coincident_hits:

        #Initialise empty array to store clusters.

        coincident_event_clustered = [coincidence_event[0],coincidence_event[1],[]]

        # coincident_event_clustered = [['Event x',TIMEBIN, [[[RPC1_PHI_CLUSTERS],[RPC1_ETA_CLUSTERS]],[[...],[...]]]]]
        
        # coincidence_event : ['Event x', TIMEBIN, [hit_locations]]
        # TIMEBIN is the initial time of the temporally clustered events.
        # hit_locations = [[RPC,CHANNEL,HIT_TIME,eta/phi]...]

        hit_locations = coincidence_event[2]
        
        #Extract hit_locations in phi and eta directions.
        phi_locations = [x for x in hit_locations if x[3]=='phi']
        eta_locations = [x for x in hit_locations if x[3]=='eta']

        #Sort by channels
        phi_locations = sorted(phi_locations, key=lambda x: x[1])
        eta_locations = sorted(eta_locations, key=lambda x: x[1])

        for RPC in range(6):

            rpc_phi_clusters = []
            rpc_eta_clusters = [] 

            i = 0
            for index,hit in enumerate([x for x in phi_locations if x[0]==RPC]):
                if index==0:
                    previous_element = hit[1]
                    rpc_phi_clusters.append([hit])
                else:
                    if abs(hit[1] - previous_element) > 1:
                        # Hit is not part of the same cluster, intiate a new cluster
                        rpc_phi_clusters.append([hit])
                        i += 1
                    else:
                        # Hit is part of the same cluster
                        rpc_phi_clusters[i].append(hit)
                    previous_element = hit[1]

            j = 0
            for index,hit in enumerate([x for x in eta_locations if x[0]==RPC]):
            
                if index == 0:
                    previous_element = hit[1]
                    rpc_eta_clusters.append([hit])
                else:
                    if abs(hit[1] - previous_element) > 1:
                        # Hit is not part of the same cluster
                        rpc_eta_clusters.append([hit])
                        j += 1
                    else:
                        # Hit is part of the same cluster
                        rpc_eta_clusters[j].append(hit)
                    previous_element = hit[1]

            rpc_combined = [rpc_phi_clusters,rpc_eta_clusters]

            coincident_event_clustered[2].append(rpc_combined)

        coincident_hits_clustered.append(coincident_event_clustered)

    return coincident_hits_clustered

def check_event_attributes(event,min_chamber_number,min_RPC_number):
    #Used in filter_events() function to decide whether or not to save an event as described by the user's inputs.

    # event = ['Event x',TIMEBIN, [[[RPC1_PHI_CLUSTERS],[RPC1_ETA_CLUSTERS]],[[...],[...]],...]

    RPC_counter = 0
    chamber_counter = 0
    condition_1 = False
    condition_2 = False
    condition_3 = False

    for RPC in range(6):
        if RPC<3:
            #Checking triplet layer.
            if event[2][RPC][0] and event[2][RPC][1]:
                #Reqiure       phi^              and               eta^       strips to go off
                RPC_counter+=1 
                #If RPC has two eta and phi strips going off then consider it "hit"
                if not condition_1:
                    #Count triplet chamber being hit.
                    chamber_counter+=1
                    condition_1 = True
        elif RPC == 3:
            #Singlet layer
            if event[2][RPC][0] and event[2][RPC][1]:
                RPC_counter+=1
                if not condition_2:
                    chamber_counter+=1
                    condition_2 = True
        else:
            #Doublet layer
            if event[2][RPC][0] and event[2][RPC][1]:
                RPC_counter+=1
                if not condition_3:
                    chamber_counter+=1
                    condition_3 = True

    return RPC_counter >= min_RPC_number and chamber_counter >= min_chamber_number

def filter_events(events,min_chamber_number,min_RPC_number):
    #Initiliase array of filtered events
    filtered_events = []

    for event in events:
        if check_event_attributes(event,min_chamber_number,min_RPC_number):
            filtered_events.append(event)

    print(f"Number of events in filter = {len(filtered_events)}")
    
    return filtered_events

#RECONSTRUCTION

def extract_coords(event,max_cluster_size):

    #This function converts spatially clusters in RPCs into x and y coordinates (z given by RPC number)
    # event = ['Event x',TIMEBIN, [[[RPC1_PHI_CLUSTERS],[RPC1_ETA_CLUSTERS]],[[...],[...]],...]

    #Extract x and y coords of cluster in event
    RPC_heights = [0.6,1.8,3.0,61.8,121.8,123] #Heights of middle point of each RPC, measured from the bottom of the Triplet Low RPC. Units are cm.
    distance_per_phi_channel = 2.7625 #cm
    distance_per_eta_channel = 2.9844 #cm
    
    coords = []


    for RPC in range(6):
        
        x_clusters = [x for x in event[2][RPC][0] if len(x)<=max_cluster_size] #phi direction
        y_clusters = [y for y in event[2][RPC][1] if len(y)<=max_cluster_size] #eta direction

        #Finding size of largest cluster, consider coordinates bad if largest cluster is larger than 6.
        x_clusters_lengths = [len(x) for x in event[2][RPC][0]]
        y_clusters_lengths = [len(y) for y in event[2][RPC][1]]

        max_length = max(max(x_clusters_lengths, default=0), max(y_clusters_lengths, default=0))

        x_coords = []
        y_coords = []

        for x_cluster in x_clusters:
           #x_cluster = [[RPC,CHANNEL,TIME,'phi'],...]
            phi_channels = [x[1] for x in x_cluster]
            #Convert the channel number into a measurement along the RPC.
            x_values = [(phi_channel+0.5)*distance_per_phi_channel for phi_channel in phi_channels]
            x_coords.append(np.mean(x_values))
            #x_error

        for y_cluster in y_clusters:
            #y_cluster = [[RPC,CHANNEL,TIME,'eta'],...]
            eta_channels_corrected = [31-y[1] for y in y_cluster] #corrected for labelling from 0 to 31.
            y_values = [(channel_num+0.5)*distance_per_eta_channel for channel_num in eta_channels_corrected]
            y_coords.append(np.mean(y_values))
            #x_error

        if x_coords and y_coords and max_length<6:
            coords.append([x_coords, y_coords])
        else:
            coords.append([[],[]])

    #coords = [[RPC1_coords],[RPC2_coords],[RPC3_coords],...]
    return(coords)

def generate_hit_coords_combo(coords, RPC_heights, combinations=None, hit_coords=None, depth=0):

    if combinations is None:
        combinations = []
    if hit_coords is None:
        hit_coords = []

    if depth == len(coords):
        combinations.append(hit_coords.copy())
        return combinations

    x_values = coords[depth][0]
    y_values = coords[depth][1]

    if not x_values or not y_values:
        return generate_hit_coords_combo(coords, RPC_heights, combinations, hit_coords, depth + 1)

    for x in x_values:
        for y in y_values:
            if x is not None and y is not None and isinstance(x, (int, float)) and isinstance(y, (int, float)):
                hit_coords.append([x, y, RPC_heights[depth]])
                generate_hit_coords_combo(coords, RPC_heights, combinations, hit_coords, depth + 1)
                hit_coords.pop()

    return combinations

def fit_event(coordinates):
    #Coordinates = [[x0,y0,z0],[x1,y1,z1],...,[x5,y5,z2]]
    #Z coordinate given by height of relevant RPC.
    #Using SVD

    centroid = np.mean(coordinates, axis=0)
    subtracted = coordinates-centroid

    # performing SVD
    _, _, V = np.linalg.svd(subtracted)
    
    # find the direction vector (which is the right singular vector corresponding to the largest singular value)
    direction = V[0, :]

    # A line is defined by the average and its direction
    p0 = centroid
    d = direction

    #Work out residuals. Minimise this to find best fit (from possible combos)

    residuals = 0

    for point in coordinates:

        x_p0 = np.subtract(point,p0)

        min_distance = np.linalg.norm(np.cross(x_p0,d)) / np.linalg.norm(d)

        residuals += (min_distance)**2

    return p0, d, residuals

def reconstruct(event,max_cluster_size):
    max_residual = 100

    # event = ['Event x',TIMEBIN, [[[RPC1_PHI_CLUSTERS],[RPC1_ETA_CLUSTERS]],[[...],[...]],...]
    RPC_heights = [0.6,1.8,3.0,61.8,121.8,123] #Heights of middle point of each RPC, measured from the bottom of the Triplet Low RPC. Units are cm.

    #Extract x and y coords of cluster in event

    coords = extract_coords(event,max_cluster_size)

     # Count the number of empty RPCs
    empty_RPC_count = sum(1 for item in coords if item == [[], []])

    # # Count the number of non-empty RPCs at the beginning of the coords list
    # non_empty_count = sum(1 for item in coords[:3] if item != [[], []])

    # # If there are only three non-empty RPCs at the beginning, exit the function
    # if empty_RPC_count >= 3 and non_empty_count == 3:
    #     print("Failed to reconstruct, not enough chambers crossed")
    #     return None  # Exit the function

    # If less than 3 elements of coords are occupied, exit the function
    if empty_RPC_count > 3:
        print("Failed to reconstruct, not enough coords")
        return None  # Exit the function

    #ITERATING OVER EVERY POSSIBLE COMBINATION OF x,y,z over all 3 RPCs (limited to one x,y per RPC).
    #Doesn't look particularly nice, but there are not many coordinates to loop over usually....

    combinations = ANT.generate_hit_coords_combo(coords,RPC_heights)

    #Now for each combo in combinations, attempt to reconstruct a path. See which one gives the best trajectory.

    #If success, print parameters of fitting function.
    #If fail, print reconstruction failed.

    residuals_current = np.inf
    optimised_coords = None
    optimised_d= None
    optimised_centroid= None

    for ind,combo in enumerate(combinations):

        centroid, d, residuals = ANT.fit_event(combo)
        if residuals < residuals_current:

            # If new fit is better than old then replace old fit.
            residuals_current = residuals
            optimised_centroid = centroid
            optimised_d = d
            optimised_coords = combinations[ind]

    #FLIP VECTOR SO IT POINTS DOWN!
    #Just minor fix for now. In future, find average time of cluster in highest plate and compare with lowest hit plate to
    #determine the vertical direction of a track.
    
    if optimised_d[2]  > 0:
        optimised_d = np.multiply(optimised_d,-1)

    if residuals_current<max_residual:
        return optimised_centroid, optimised_d, optimised_coords, combinations, residuals_current
    else:
        print("Failed to reconstruct")
        #return optimised_centroid, optimised_d, optimised_coords, combinations, residuals_current
        return None

def line_equation(t, a, b, c, x0, y0, z0):
    # Define the parametric equation of the line
    x = x0 + a*t
    y = y0 + b*t
    z = z0 + c*t
    return x, y, z

def non_interactive_muon_plot(centroid,d,event_coords):
    #Coefficients = [a,b,c]
    #event_coords = [[x0,y0,z0],[x1,y1,z1],...,[x5,y5,z5]]

   # Generate line coordinates
    t_values = np.linspace(-100, 100, 100)
    line_coordinates = centroid.reshape(-1, 1) + d.reshape(-1, 1) * t_values.reshape(1, -1)

    x_coords = [x[0] for x in event_coords]
    y_coords = [y[1] for y in event_coords]
    z_coords = [z[2] for z in event_coords]


    # Create a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the data points
    ax.scatter(x_coords, y_coords, z_coords, color='r', marker='o', label='Data Points')

    # Plot the line
    ax.plot(*line_coordinates, color='b', label='Fitted Line')

    # Plot centroid
    ax.scatter(centroid[0], centroid[1], centroid[2], color='green', label='Centroid')

    # Set labels
    ax.set_xlabel('X/cm')
    ax.set_ylabel('Y/cm')
    ax.set_zlabel('Z/cm')

    plt.xlim(0,180)
    plt.ylim(0,99)
    ax.set_zlim(-10,130)

    ax.view_init(-140, 60)

    # Add a legend
    ax.legend()

    plt.show()

def interactive_muon_plot(centroid,d,event_coords):
    #Coefficients = [a,b,c]
    #event_coords = [[x0,y0,z0],[x1,y1,z1],...,[x5,y5,z5]]

    # Generate line coordinates
    t_values = np.linspace(-150, 150, 100)
    line_coordinates = centroid.reshape(-1, 1) + d.reshape(-1, 1) * t_values.reshape(1, -1)

    x_coords = [x[0] for x in event_coords]
    y_coords = [y[1] for y in event_coords]
    z_coords = [z[2] for z in event_coords]


    RPC_origins = [[0,0,0],[0,0,1.2],[0,0,2.4],[0,0,61.2],[0,0,121.2],[0,0,122.4]]
    RPC_dimensions = [180,99,1.2]

    #Calculate vertices for RPCs
    rpc_vertices = []
    for i in range(6):
        rpc_vertices.append(calculate_cuboid_vertices(RPC_origins[i],RPC_dimensions))

    # Configure Plotly to be rendered inline in the notebook.
    plotly.offline.init_notebook_mode()

    # Configure the trace for the RPCs
    rpc_0 = go.Mesh3d(
        x=[vertex[0] for vertex in rpc_vertices[0]],
        y=[vertex[1] for vertex in rpc_vertices[0]],
        z=[vertex[2] for vertex in rpc_vertices[0]],
        i=[0, 1, 2, 3, 0, 4, 5, 6, 7, 4, 5, 1, 2, 6, 7, 3],  # Indices for faces
        j=[1, 2, 3, 0, 4, 5, 6, 7, 5, 6, 2, 3, 7, 4, 0, 1],  # Indices for faces
        k=[2, 3, 0, 1, 6, 7, 4, 5, 1, 5, 6, 2, 3, 7, 4, 0],  # Indices for faces
        opacity=0.1,  # Set the opacity as needed
        color='green'  # Set the color of the cuboid
    )

    rpc_1 = go.Mesh3d(
        x=[vertex[0] for vertex in rpc_vertices[1]],
        y=[vertex[1] for vertex in rpc_vertices[1]],
        z=[vertex[2] for vertex in rpc_vertices[1]],
        i=[0, 1, 2, 3, 0, 4, 5, 6, 7, 4, 5, 1, 2, 6, 7, 3],  # Indices for faces
        j=[1, 2, 3, 0, 4, 5, 6, 7, 5, 6, 2, 3, 7, 4, 0, 1],  # Indices for faces
        k=[2, 3, 0, 1, 6, 7, 4, 5, 1, 5, 6, 2, 3, 7, 4, 0],  # Indices for faces
        opacity=0.1,  # Set the opacity as needed
        color='green'  # Set the color of the cuboid
    )

    rpc_2 = go.Mesh3d(
        x=[vertex[0] for vertex in rpc_vertices[2]],
        y=[vertex[1] for vertex in rpc_vertices[2]],
        z=[vertex[2] for vertex in rpc_vertices[2]],
        i=[0, 1, 2, 3, 0, 4, 5, 6, 7, 4, 5, 1, 2, 6, 7, 3],  # Indices for faces
        j=[1, 2, 3, 0, 4, 5, 6, 7, 5, 6, 2, 3, 7, 4, 0, 1],  # Indices for faces
        k=[2, 3, 0, 1, 6, 7, 4, 5, 1, 5, 6, 2, 3, 7, 4, 0],  # Indices for faces
        opacity=0.1,  # Set the opacity as needed
        color='green'  # Set the color of the cuboid
    )

    rpc_3 = go.Mesh3d(
        x=[vertex[0] for vertex in rpc_vertices[3]],
        y=[vertex[1] for vertex in rpc_vertices[3]],
        z=[vertex[2] for vertex in rpc_vertices[3]],
        i=[0, 1, 2, 3, 0, 4, 5, 6, 7, 4, 5, 1, 2, 6, 7, 3],  # Indices for faces
        j=[1, 2, 3, 0, 4, 5, 6, 7, 5, 6, 2, 3, 7, 4, 0, 1],  # Indices for faces
        k=[2, 3, 0, 1, 6, 7, 4, 5, 1, 5, 6, 2, 3, 7, 4, 0],  # Indices for faces
        opacity=0.1,  # Set the opacity as needed
        color='green'  # Set the color of the cuboid
    )

    rpc_4 = go.Mesh3d(
        x=[vertex[0] for vertex in rpc_vertices[4]],
        y=[vertex[1] for vertex in rpc_vertices[4]],
        z=[vertex[2] for vertex in rpc_vertices[4]],
        i=[0, 1, 2, 3, 0, 4, 5, 6, 7, 4, 5, 1, 2, 6, 7, 3],  # Indices for faces
        j=[1, 2, 3, 0, 4, 5, 6, 7, 5, 6, 2, 3, 7, 4, 0, 1],  # Indices for faces
        k=[2, 3, 0, 1, 6, 7, 4, 5, 1, 5, 6, 2, 3, 7, 4, 0],  # Indices for faces
        opacity=0.1,  # Set the opacity as needed
        color='green'  # Set the color of the cuboid
    )

    rpc_5 = go.Mesh3d(
        x=[vertex[0] for vertex in rpc_vertices[5]],
        y=[vertex[1] for vertex in rpc_vertices[5]],
        z=[vertex[2] for vertex in rpc_vertices[5]],
        i=[0, 1, 2, 3, 0, 4, 5, 6, 7, 4, 5, 1, 2, 6, 7, 3],  # Indices for faces
        j=[1, 2, 3, 0, 4, 5, 6, 7, 5, 6, 2, 3, 7, 4, 0, 1],  # Indices for faces
        k=[2, 3, 0, 1, 6, 7, 4, 5, 1, 5, 6, 2, 3, 7, 4, 0],  # Indices for faces
        opacity=0.1,  # Set the opacity as needed
        color='green'  # Set the color of the cuboid
    )

    # Configure the trace.
    trace = go.Scatter3d(
        x=x_coords,  # <-- Put your data instead
        y=y_coords,  # <-- Put your data instead
        z=z_coords, # <-- Put your data instead
        mode='markers',
        marker={
            'size': 10,
            'opacity': 1.0,
        }
    )

    #Plot zenith angle
    # zenith = go.Vector(...)

    # Extract x, y, z coordinates from line coordinates
    x_line = line_coordinates[0]
    y_line = line_coordinates[1]
    z_line = line_coordinates[2]
    
    # Configure the trace for the line
    trace_line = go.Scatter3d(
        x=x_line,
        y=y_line,
        z=z_line,
        mode='lines',
        line={
            'color': 'red',
            'width': 2,
        },
        name='Line'
    )

    trace_centroid = go.Scatter3d(
        x=[centroid[0]],
        y=[centroid[1]],
        z=[centroid[2]],
        mode='markers',
        marker={
            'size': 5,
            'color': 'green',
        },
        name='Centroid'
    )

    # Configure the layout.
    layout = go.Layout(
        margin={'l': 0, 'r': 0, 'b': 0, 't': 0},
        scene=dict(
            xaxis=dict(showgrid=True, gridcolor='rgb(211,211,211)', gridwidth=2,range=[-5,185], title = "x/cm"),
            yaxis=dict(showgrid=True, gridcolor='rgb(211,211,211)', gridwidth=2,range=[-5, 105],title = "y/cm"),
            zaxis=dict(showgrid=True, gridcolor='rgb(211,211,211)', gridwidth=2,range=[-5, 125],title="z/cm")
        )
    )

    #Include trace_centroid in data= [] to plot centroid on the plot.

    data = [trace,trace_line,rpc_0,rpc_1,rpc_2,rpc_3,rpc_4,rpc_5]

    plot_figure = go.Figure(data=data, layout=layout)

    # Render the plot.
    plotly.offline.iplot(plot_figure)

def calculate_cuboid_vertices(origin, dimensions):
    # Function to calculate vertices of cuboid from origin and dimensions
    x_min, y_min, z_min = origin
    x_max = x_min + dimensions[0]
    y_max = y_min + dimensions[1]
    z_max = z_min + dimensions[2]

    vertices = [
        [x_min, y_min, z_min], [x_max, y_min, z_min],
        [x_max, y_max, z_min], [x_min, y_max, z_min],
        [x_min, y_min, z_max], [x_max, y_min, z_max],
        [x_max, y_max, z_max], [x_min, y_max, z_max]
    ]

    return vertices

#PLOTTING FUNCTIONS

def plot_dark_clustering(phi_cluster_distribution,eta_cluster_distribution, time ="60 Seconds"):
    
    RPC_description = ['Triplet Low','Triplet Mid','Triplet Top','Singlet','Doublet Low','Doublet Top']
    multi_array = []
    current_directory=  os.path.dirname(os.getcwd())

    for RPC in range(len(RPC_description)):
        if phi_cluster_distribution[RPC]:

            fig, ax = plt.subplots(1, figsize=(16, 8), dpi=100)
            bin_edges = [i - 0.5 for i in range(min(phi_cluster_distribution[RPC]), max(phi_cluster_distribution[RPC]) + 2)]
            plt.hist(phi_cluster_distribution[RPC], bins=bin_edges, align='mid', rwidth=0.5)
            plt.xlabel('Cluster size')
            plt.ylabel('Count')
            plt.xlim(0, max(phi_cluster_distribution[RPC]) + 5)
            plt.title(f'Cluster distributions, Phi Direction, RPC {RPC_description[RPC]}')
            plt.tight_layout()
            plt.savefig(current_directory+f"\\Figures\\Phi{RPC}.png")
            multi_array.append(current_directory+f"\\Figures\\Phi{RPC}.png")

    for RPC in range(len(RPC_description)):
        if eta_cluster_distribution[RPC]:

            fig, ax = plt.subplots(1, figsize=(16, 8), dpi=100)
            bin_edges = [i - 0.5 for i in range(min(eta_cluster_distribution[RPC]), max(eta_cluster_distribution[RPC]) + 2)]
            plt.hist(eta_cluster_distribution[RPC], bins=bin_edges, align='mid', rwidth=0.5)
            plt.xlabel('Cluster size')
            plt.ylabel('Count')
            plt.xlim(0, max(eta_cluster_distribution[RPC]) + 5)
            plt.title(f'Cluster distributions, Eta Direction, RPC {RPC_description[RPC]}')
            plt.tight_layout()  # Adjust layout
            plt.savefig(current_directory+f"\\Figures\\Eta{RPC}.png")
            multi_array.append(current_directory+f"\\Figures\\Eta{RPC}.png")

    anPlot.combinePlots(multi_array,f"Dark Clustering, {time}")
    for plot in multi_array:
        os.remove(plot)

def plot_event_cluster_COMBINED(event_cluster):
    #PLOTS ALL STRIPS FIRED IN A EVENT_CLUSTER ON THE SAME RPC!

    #event_cluster = ['Event X', cluster_initial_time, [hit_locations]]
    #hit_locations = [[RPC,Channel,Hit_time,'eta/phi'],...]

    initial_time = event_cluster[1]
    RPC_events = event_cluster[2]

    # Extracting x and y data

    # [x[1],x[2]] = [Channel, Hit_time]
    phi_data = [[x[1],x[2]] for x in RPC_events if x[3]=='phi']
    eta_data = [[y[1],y[2]] for y in RPC_events if y[3]=='eta']

    # Creating histogram
    plt.figure(figsize=(12, 8))
    # Plotting bars for active phi strips
    for x in phi_data:
        plt.bar(x[0], 32, width=0.8, color='green', alpha=0.5)
        plt.text(x[0]-0.5,16,f"{x[1]} ns",rotation=90)

    for y in eta_data:
        plt.barh(31-y[0], 64, height=0.8, color='red', alpha=0.5)
        plt.text(32,(31-y[0])-0.25,f"{x[1]} ns")

    #plt.colorbar(label='Counts')
    plt.xlabel('Phi Channel')
    plt.ylabel('Eta Channel')
    plt.title(f'Fired channels for {event_cluster[0]} at {initial_time} ns COMBINED RPC PLOT')
    plt.tight_layout()

    plt.xlim(-1, 64)
    plt.ylim(-1, 32)

    # Set x-axis ticks every 2 units
    ax = plt.gca()
    ax.xaxis.set_major_locator(MultipleLocator(2))

    # Set y-axis ticks to cover the entire range of eta channels
    ax.set_yticks(range(32))

    # Reverse y-axis ticks to match eta channel numbers
    ax.set_yticklabels(range(31, -1, -1))

    
    #plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_event_cluster(event_cluster):
    #PLOT HIT MAP FOR EVENT_CLUSTER IN AN EVENT

    RPC_description = ['Triplet Low','Triplet Mid','Triplet Top','Singlet','Doublet Low','Doublet Top']

    #event_cluster = ['Event X', cluster_initial_time, [hit_locations]]
    #hit_locations = [[RPC,Channel,Hit_time,'eta/phi'],...]

    initial_time = event_cluster[1]
    RPC_events = event_cluster[2]

    # Extracting x and y data

    # [x[1],x[2]] = [Channel, Hit_time]
    phi_data = [[x[1],x[2]] for x in RPC_events if x[3]=='phi']
    eta_data = [[y[1],y[2]] for y in RPC_events if y[3]=='eta']

    for RPC in range(6):

        phi_data = [[x[1],x[2]] for x in RPC_events if x[3]=='phi' and x[0]==RPC]
        eta_data = [[y[1],y[2]] for y in RPC_events if y[3]=='eta' and y[0]==RPC]

        # Creating histogram
        plt.figure(figsize=(12, 8))
        # Plotting bars for active phi strips
        for x in phi_data:
            plt.bar(x[0], 32, width=0.8, color='green', alpha=0.5)
            plt.text(x[0]-0.5,16,f"{x[1]} ns",rotation=90)

        for y in eta_data:
            plt.barh(31-y[0], 64, height=0.8, color='red', alpha=0.5)
            plt.text(32,(31-y[0])-0.25,f"{y[1]} ns")

        #plt.colorbar(label='Counts')
        plt.xlabel('Phi Channel')
        plt.ylabel('Eta Channel')
        plt.title(f'Fired channels for {event_cluster[0]} at {initial_time} ns, RPC {RPC_description[RPC]}')
        plt.tight_layout()

        plt.xlim(-1, 64)
        plt.ylim(-1, 32)

        # Set x-axis ticks every 2 units
        ax = plt.gca()
        ax.xaxis.set_major_locator(MultipleLocator(2))

        # Set y-axis ticks to cover the entire range of eta channels
        ax.set_yticks(range(32))

        # Reverse y-axis ticks to match eta channel numbers
        ax.set_yticklabels(range(31, -1, -1))

        
        #plt.grid(True)
        plt.tight_layout()
        plt.show()

def plot_event_cluster_save(event_cluster):

    multi_array = []
    current_directory=  os.path.dirname(os.getcwd())

    RPC_description = ['Triplet Low','Triplet Mid','Triplet Top','Singlet','Doublet Low','Doublet Top']

    #event_cluster = ['Event X', cluster_initial_time, [hit_locations]]
    #hit_locations = [[RPC,Channel,Hit_time,'eta/phi'],...]

    initial_time = event_cluster[1]
    RPC_events = event_cluster[2]

    # Extracting x and y data

    # [x[1],x[2]] = [Channel, Hit_time]
    phi_data = [[x[1],x[2]] for x in RPC_events if x[3]=='phi']
    eta_data = [[y[1],y[2]] for y in RPC_events if y[3]=='eta']

    for RPC in range(6):

        phi_data = [[x[1],x[2]] for x in RPC_events if x[3]=='phi' and x[0]==RPC]
        eta_data = [[y[1],y[2]] for y in RPC_events if y[3]=='eta' and y[0]==RPC]

        # Creating histogram
        plt.figure(figsize=(12, 8))
        # Plotting bars for active phi strips
        for x in phi_data:
            plt.bar(x[0], 32, width=0.8, color='green', alpha=0.5)
            plt.text(x[0]-0.5,16,f"{x[1]} ns",rotation=90)

        for y in eta_data:
            plt.barh(31-y[0], 64, height=0.8, color='red', alpha=0.5)
            plt.text(32,(31-y[0])-0.25,f"{y[1]} ns")

        #plt.colorbar(label='Counts')
        plt.xlabel('Phi Channel')
        plt.ylabel('Eta Channel')
        plt.title(f'Fired channels for {event_cluster[0]} at {initial_time} ns, RPC {RPC_description[RPC]}')
        plt.tight_layout()

        plt.xlim(-1, 64)
        plt.ylim(-1, 32)

        # Set x-axis ticks every 2 units
        ax = plt.gca()
        ax.xaxis.set_major_locator(MultipleLocator(2))

        # Set y-axis ticks to cover the entire range of eta channels
        ax.set_yticks(range(32))

        # Reverse y-axis ticks to match eta channel numbers
        ax.set_yticklabels(range(31, -1, -1))

        #plt.grid(True)
        plt.tight_layout()
        plt.savefig(current_directory+f"\\Figures\\{event_cluster[0]}{RPC}.png")
        multi_array.append(current_directory+f"\\Figures\\{event_cluster[0]}{RPC}.png")

    anPlot.combinePlots(multi_array,f"Hit_map,{event_cluster[0]},{initial_time}ns")
    for plot in multi_array:
        os.remove(plot)

def combinePlots(plots,imname):
    images = [Image.open(x) for x in plots]
    widths, heights = zip(*(i.size for i in images))

    total_width = int(2*widths[0])
    if(len(plots)|2>0):
        max_height = int((sum(heights)+heights[0])/2)
    else:
        max_height = int(sum(heights)/2)

    new_im = Image.new('RGB', (total_width, max_height))

    x_offset = 0
    y_offset = 0
    even = True
    for im in images:
        if even:
            new_im.paste(im, (x_offset,y_offset))
            x_offset += im.size[0]
            even = False
        else:
            new_im.paste(im,(x_offset,y_offset))
            x_offset = 0
            y_offset += im.size[1]
            even = True

    new_im.save(imname.strip(" ")+'.pdf')

def heatFromFile(dataFile, time=240, name="HeatMap"):
    #Plots heat maps from triggered data, showing the hit rate in each rpc channel. 2D plots designed to replicate RPC layout and channel counting direction.
    thisData = importDatafile(dataFile)
    thisHitData = {}
    addresses = ['ee00','ee01','ee02','ee03','ee04']
    for tdc in range(5):
        thisHitData[addresses[tdc]] = countChannels(thisData[tdc])
    anPlot.makeHitMaps(thisHitData,name,False,unit='hz',time=time)

def plotHitCounts(histogram, name):
    #Plot the number of hits per event to determine the mean hits within any given TDC, as well as see the long tails from correlated noise
    fig, ax = plt.subplots(1, figsize=(6, 4), dpi=100)
    lab = hep.atlas.label(com=False,data=True, label="Internal")
    lab[2].set_text(" ")
    hep.histplot(np.histogram(histogram,bins=[x-0.5 for x in range(140)]), label='TDC 0')
    plt.xlabel('TDC Hits')
    plt.ylabel('Events')
    plt.title(name)
    plt.xlim([-0.4,40.5])
    #plt.ylim([1,3000])
    plt.yscale('log')
    plt.savefig(name.strip(" ")+"chanCounts.png")
    return name.strip(" ")+"chanCounts.png"

def plotEventTimes(inputData, name):
    fig, ax = plt.subplots(1, figsize=(6, 4), dpi=100)
    lab = hep.atlas.label(com=False,data=True, label="Internal")
    lab[2].set_text(" ")
    hep.histplot(np.histogram(inputData,bins=[x-0.5 for x in range(1000)]), label='TDC 0')
    plt.xlabel('TDC Hit Times')
    plt.ylabel('Hits')
    plt.title(name)
    #plt.xlim([-0.4,40.5])
    #plt.ylim([1,3000])
    plt.yscale('log')
    plt.savefig(name.strip(" ")+"hitTimes.png")
    return name.strip(" ")+"hitTimes.png"

def convert_cluster_to_plot(clustered_event):

     #Just convert the clustered_event into the old event_cluster format to be able to plot it using old fucntion.

     # clustered_event = ['Event x',TIMEBIN, [[[RPC1_PHI_CLUSTERS],[RPC1_ETA_CLUSTERS]],[...]]]
    
    event_cluster = [clustered_event[0],clustered_event[1]]

    hit_locations = []

    for RPC in clustered_event[2]:
        for phi_cluster in RPC[0]:
            for hit in phi_cluster:
                hit_locations.append(hit)
        for eta_cluster in RPC[1]:
            for hit in eta_cluster:
                hit_locations.append(hit)

    event_cluster.append(hit_locations)

    return event_cluster

def plot_zenith_angle_distribution(angles, title):

    #Histogram plot given anuglar information of hits. angles is just array with elements that are angles from 0 to pi/2 radians.
    num_points = len(angles)

    theta_vals = np.linspace(0, np.pi / 2, 100000, endpoint=False)
    probs = [np.sin(x) * (np.cos(x))**2 for x in theta_vals]
    norm_probs = np.multiply(1 / (np.sum(probs)), probs)

    cdf = np.cumsum(norm_probs)
    cdf_spacing = np.pi / 2 / 1e5

    plt.figure(figsize=(16,8))

    # Plot histogram of data.
    num_bins = 100  # Adjust the number of bins as needed
    counts, bin_edges, _ = plt.hist(angles, bins=num_bins, alpha=0.7, edgecolor='black', label='ProANUBIS muon distribution')

    bin_midpoints = [(bin_edges[i] + bin_edges[i + 1]) / 2 for i in range(len(bin_edges) - 1)]
    midpoint_freq = []

    for j in range(len(bin_edges) - 1):
        start_bin_index = int(np.floor(bin_edges[j] / cdf_spacing))
        end_bin_index = int(np.floor(bin_edges[j + 1] / cdf_spacing))
        cum_prob = cdf[min(end_bin_index, len(cdf) - 1)] - cdf[start_bin_index]
        freq = cum_prob * num_points
        midpoint_freq.append(freq)

    # Plot distribution from literature.
    plt.plot(bin_midpoints, midpoint_freq, c='red', label='Cosmic Ray muon distribution')

    plt.xlabel("Zenith Angle/ degrees")  # Change the x-axis label
    plt.ylabel("Count")

    plt.legend(loc='upper right')
    
    # Adjust x-axis limits
    plt.xlim(left=0, right=np.pi / 2)
  
    # Convert x-axis ticks from radians to degree
    x_ticks_degrees = np.linspace(0, 90, num=10)
    x_ticks_radians = np.radians(x_ticks_degrees)

    plt.xticks(x_ticks_radians, x_ticks_degrees)

    plt.title(f"{title}")

    plt.show()

def plot_angle_distribution(angles, title):
    #Plot angular distribution given phi and eta angular distribution separately.

    plt.figure(figsize=(16,10))

    # Plot histogram with counts
    plt.hist(angles, bins=61, density=True,edgecolor='black',alpha=0.7, label='ProANUBIS muon distribution')

    # Define the cosine squared function and plot it
    x = np.linspace(-np.pi/2, np.pi/2, 100)
    cos_squared = 2/np.pi * np.cos(x)**2
    plt.plot(x, cos_squared, color='red', label='Cosmic Muon distribution')

    # Convert radians to degrees for x-ticks
    x_ticks_degrees = np.linspace(-90, 90, num=19)
    x_ticks_radians = np.radians(x_ticks_degrees)

    # Set x-ticks labels and positions
    plt.xticks(x_ticks_radians, x_ticks_degrees)

    plt.text(-1,1.1,f"Number of events = {len(angles)}")

    # Customize the plot
    plt.xlabel('Angle/ degrees from axis perpendicular to surface of RPCs')
    plt.ylabel('Relative Occurence')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.show()

def extract_angles_phi_eta(filtered_events):

    #Input is filtered_events, output of ANT.filter_events() function

    angles_eta = []
    angles_phi = []

    for i,filtered_event in enumerate(filtered_events):

        result = reconstruct(filtered_event,3)

        if result is not None:
            #Only save angles that actually were reconstructed well
            
            # a.b = |a||b|cos(x)

            #eta angle. 
            #work out the projection of the direction vector in the plane.
            
            v_parr_eta = np.array([0,result[1][1],result[1][2]])

            theta_eta = np.arccos(np.dot(v_parr_eta,[0,0,1]) / np.linalg.norm(v_parr_eta))

            if theta_eta > np.pi / 2:
                theta_eta= np.pi - theta_eta
            
            if v_parr_eta[1] > 0:
                theta_eta*=-1

            angles_eta.append(theta_eta)

            # Phi angles
            #work out the projection of the direction vector in the plane.
            
            v_parr_phi = np.array([result[1][0],0,result[1][2]])

            theta_phi = np.arccos(np.dot(v_parr_phi,[0,0,1]) / np.linalg.norm(v_parr_phi))

            if theta_phi > np.pi / 2:
                theta_phi= np.pi - theta_phi
            
            if v_parr_phi[0] < 0:
                theta_phi*=-1

            angles_phi.append(theta_phi)

    return angles_eta, angles_phi

        #ProAnubis setup is at 45 degrees to vertical. 
        #Project direction vector onto planes to work out phi and eta angular distributions. Should be no assymmetry in phi.
        #Expect asymmetry in eta. 

