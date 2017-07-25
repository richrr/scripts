#http://stackoverflow.com/questions/11098189/awk-extract-columns-from-file-based-on-header-selected-from-2nd-file

DATAFILE=${1:-data.txt}
COLUMNFILE=${2:-list.txt}

awk -F "\t" -v colsFile="$COLUMNFILE" '
   BEGIN {
     j=1
     while ((getline < colsFile) > 0) {
        col[j++] = $1
     }
     n=j-1;
     close(colsFile)
     for (i=1; i<=n; i++) s[col[i]]=i
   }
   NR==1 {
     for (f=1; f<=NF; f++)
       if ($f in s) c[s[$f]]=f
     next
   }
   { sep=""
     for (f=1; f<=n; f++) {
       printf("%c%s",sep,$c[f])
       sep=FS
     }
     print ""
   }
' "$DATAFILE"
