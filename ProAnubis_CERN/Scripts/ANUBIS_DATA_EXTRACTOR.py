import ROOT

def importFromROOTFile(filename):
    inputRoot = ROOT.TFile(filename)
    eventTree = inputRoot.Get("Events")    
    data = [[] for tdc in range(5)]
    for event in range(eventTree.GetEntries()):
        eventTree.GetEvent(event)
        tdc = eventTree.TDC
        thisEvent = []
        hits = eventTree.Words
        for hit in hits:
            thisEvent.append(hit)
        data[tdc].append(thisEvent)
        #eventTime = eventTree.TimeStamp.AsString()
    endTime = inputRoot.Get("endTime").AsString()
    return data, endTime

file_directory = r"File directory"

data, endTime = importFromROOTFile(file_directory)

DataFile = open('/DataFile.txt', 'w')
DataFile.writelines(data)
DataFile.close()