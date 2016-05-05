push:
	-git commit -am "Cambios"
	-git push origin master

pull:
	-git pull

clean:
	find . -name "*.pyc" -and -name "*~" -exec rm {} \;
