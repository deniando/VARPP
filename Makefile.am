EXTRA_DIST = autogen.sh


SUBDIRS =  dev


count:
	line=$$(find . \( -name "*.[c,h,R]" -o -name "*.am" -o -name "*.org" -o -name "*.sh" \)  -type f -print0  | xargs -0 wc -l | grep total); echo	$$line
