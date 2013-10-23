PROGRAMS = $(shell ls *.py)
DROPBOX_FOLDER = /Users/will/Dropbox/JorgeBowshocks/LuisLL/programas

all: install_to_dropbox
install_to_dropbox: 
	rsync -av $(PROGRAMS) $(DROPBOX_FOLDER)
