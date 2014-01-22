#Convert NEWS.md to what E wants.  Shamelessly stolen from Yihui

NEWS     = NEWS

all: news

news: $(NEWS)
        sed -e 's/^-/  -/' -e 's/^## *//' -e 's/^#/\t\t/' <NEWS.md | fmt -80 >$(NEWS)