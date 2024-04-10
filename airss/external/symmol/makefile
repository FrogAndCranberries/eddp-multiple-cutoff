.DEFAULT:
	-touch $@
all: symmol
symmol: symmol.f 
	-wget -c --timeout=10 --tries=1 -nv https://www.mtg.msm.cam.ac.uk/files/symmol.zip
	-unzip -o symmol.zip symmol.f
	-patch symmol.f symmol.patch
	-$(FC) -w -o symmol symmol.f
clean:
	-rm -f VERSION_TXT example.* o.* idee.doc symmol.man symmol.f symmol
