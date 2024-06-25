import synapseclient
import synapseutils
import sys

sys.stdout = open('synapse.txt','wt')

try:

  syn = synapseclient.Synapse()

  syn.login("djuna@mit.edu", "PW")

  synapseutils.sync.syncToSynapse(syn, "manifest.txt", dryRun=False)

except Exception as e:

  print("Failed to upload!")
  print(e)
