#!/bin/bash

### these commands were used to create the BAF/LRR database

# 1. create user key
head -n1 ALL_BAF_logRR.txt | awk -F'\t' '{for(i=7; i<=NF; i+=4){split($i, n, "."); print ++u, n[1]}}' > DB_userkey.txt

# 2. create snp key
awk '{print s++, $1, $3, $4}' ALL_BAF_logRR.txt > DB_snpkey.txt
awk 'NR>1{print $1, $2, $4, $5}' DB_snpkey.txt > DB_snpkey2.txt
mv DB_snpkey2.txt DB_snpkey.txt
awk 'FNR==NR{a[$1]=$5; next} {print $1, $2, $3, $4, a[$2]}' ~/Documents/MoBa_BAF/hg18to38positions.pos DB_snpkey.txt > DB_snpkey_hg38.txt

# 3. in SQL, create tables:
CREATE TABLE ukey (uid SMALLINT NOT NULL,
	chipid CHAR(17),
	PRIMARY KEY (uid));
CREATE TABLE snpkey (snpid MEDIUMINT NOT NULL,
	rsid VARCHAR(11),
	chr TINYINT,
	pos INT,
	PRIMARY KEY (snpid));
CREATE INDEX ix_chip ON ukey(chipid);
## CREATE INDEX ix_chr ON snpkey(chr);
mysql -uroot -p moba_baf -e "LOAD DATA LOCAL INFILE 'DB_userkey.txt' INTO TABLE ukey FIELDS TERMINATED BY ' ';"
mysql -uroot -p moba_baf -e "LOAD DATA LOCAL INFILE 'DB_snpkey_hg38.txt' INTO TABLE snpkey FIELDS TERMINATED BY ' ';"
## set innodb_doublewrite=0 in mysqld.conf

## in a loop:
# 4. take a slice of 10K snps
read -s -p "db pass:" pp
for c in {11..1}
do
	echo "processing chr ${c}/22"
	awk -F'\t' -v c=${c} '
		FNR==NR{s[$2]=$1; next}
		$3==c{for(i=7; i<=NF; i+=4){
			print (i-3)/4, s[$1], $i, $(i+1)
	}}' DB_snpkey.txt ALL_BAF_logRR.txt > SLICE_BAF_logRR.txt
	echo "loading chr ${c}/22"
	mysql -uroot -p${pp} moba_baf -e \
	"SET GLOBAL innodb_buffer_pool_size=24576000000;
	SET autocommit=0; SET foreign_key_checks=0;
	CREATE TABLE raw${c} (uid SMALLINT NOT NULL,
		snpid MEDIUMINT NOT NULL,
		baf FLOAT,
		lrr FLOAT,
		CONSTRAINT snpjoiner${c} FOREIGN KEY (snpid)
		  REFERENCES snpkey(snpid),
		CONSTRAINT ujoiner${c} FOREIGN KEY (uid)
		  REFERENCES ukey(uid))
		ROW_FORMAT=COMPRESSED;
	LOAD DATA LOCAL INFILE 'SLICE_BAF_logRR.txt'
		INTO TABLE raw${c} FIELDS TERMINATED BY ' ';
	COMMIT; SET foreign_key_checks=1;"

done

## later:
ALTER TABLE raw${c} ADD PRIMARY KEY (uid, snpid);
## set innodb_doublewrite back to 1 in mysqld.conf
