# Class 16
Gonzalez A16745338

This is a simple text editor called unix Same key unix cmds so far

pwd: print working dir cd: change dir ls: list out files mkdir: make a
new dir/folder nano: a basic text editor

cp notes.txt mycopy.txt ls mv mycopy.txt somesillyname.txt ls ^q

ssh -i “~/Downloads/bimm143_gag002.pem”
ubuntu@ec2-35-165-223-67.us-west-2.compute.amazonaws.com

scp -i “~/Downloads/bimm143_gag002.pem”
ubuntu@ec2-35-165-223-67.us-west-2.compute.amazonaws.com

scp -i “~/Downloads/bimm143_gag002.pem”
ubuntu@ec2-35-165-223-67.us-west-2.compute.amazonaws.com:~/work/mm-second.x.zebrafish.tsv
myresult.tsv

``` r
c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
```

     [1] "qseqid"   "sseqid"   "pident"   "length"   "mismatch" "gapopen" 
     [7] "qstart"   "qend"     "sstart"   "send"     "evalue"   "bitscore"

``` r
## Asuming your blast results are stored in an object called 'b'
# plot(b$pident  * (b$qend - b$qstart), b$bitscore)
```

``` r
library(ggplot2)
# ggplot(b, aes(pident, bitscore)) + geom_point(alpha=0.1) 
```

``` r
# ggplot(b, aes((b$pident * (b$qend - b$qstart)), bitscore)) + geom_point(alpha=0.1) + geom_smooth()
```
