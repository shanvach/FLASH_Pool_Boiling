#!/bin/bash
for GEOFILE in geo.*
do
        echo $GEOFILE to $GEOFILE.plt
        preplot $GEOFILE $GEOFILE.plt
done
for PARTFILE in parts.*
do
        echo $PARTFILE to $PARTFILE.plt
	preplot $PARTFILE $PARTFILE.plt
done
