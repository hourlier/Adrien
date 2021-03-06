import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE(s)\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from larlite import larlite as fmwk

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
for x in xrange(len(sys.argv)-1):
    my_proc.add_input_file(sys.argv[x+1])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)
#my_proc.set_io_mode(fmwk.storage_manager.kBOTH)


# Specify output root file name
my_proc.set_ana_output_file("output_analysis trees.root");
my_proc.set_output_file("output_larlite.root");

# Attach an analysis unit ... here we use a base class which does nothing.
# Replace with your analysis unit if you wish.
my_proc.add_process(fmwk.AStarTracker())
#my_proc.set_data_to_write(larlite.event_track,"ev_track")
my_proc.get_process(0).set_producer("pandoraNu","chstatus")
my_proc.get_process(0).SetCompressionFactors(1,6) # wire,time
my_proc.get_process(0).SetVerbose(0);
my_proc.get_process(0).SetDrawOutputs(0);

print
print  "Finished configuring ana_processor. Start event loop!"
print

# Let's run it.
my_proc.run(0,2000);

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
