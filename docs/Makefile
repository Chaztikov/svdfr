DOCS=index download feedback reaction-diffusion 2dChannelFlow 2d3cChannelFlow 3dchannelflow 2dviscoelastic beam  citation algo implement
DOCS1=$(addprefix jemdocs/, $(DOCS))
HDOCS=$(addsuffix .html, $(DOCS))
#PHDOCS=$(addprefix html/, $(HDOCS))
PHDOCS=$(HDOCS)

.PHONY : docs
docs : $(PHDOCS)

.PHONY : update
update : $(PHDOCS)
	@echo -n 'Copying to server...'
	# insert code for copying to server here.
	@echo ' done.'

%.html : jemdocs/%.jemdoc MENU
	python2.7 jemdoc.py -c bibover.conf -o $@ $<

.PHONY : clean
clean :
	rm -f *.html
