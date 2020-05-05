DOCS=index beam 3dchannelflow 2dviscoelastic reaction-diffusion citation feedback algo implement
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
