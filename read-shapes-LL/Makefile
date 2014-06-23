PROGRAMS = $(shell ls *.py)
DROPBOX_FOLDER = ~/Dropbox/LuisBowshocks/programas

all: install_to_dropbox make_read_only

install_to_dropbox: 
	@rsync -av $(PROGRAMS) $(DROPBOX_FOLDER)

make_read_only: 
	@for f in $(PROGRAMS); do chmod a-w $(DROPBOX_FOLDER)/$$f; done
