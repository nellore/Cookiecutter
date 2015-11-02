#!/bin/bash

set -eu
set -o pipefail

if [ $# -lt 1 ]; then
	echo Please specify the output file name.
	echo For example: ./create_package ~/Desktop/cookiecutter_osx.tar.gz 
	exit 1
fi

OUTPUT_FILE="$1"
PACKAGE_NAME=$( basename "$OUTPUT_FILE" )
PACKAGE_NAME=$( echo "$PACKAGE_NAME" | sed "s/\..*//" )

TEMPDIR=$( mktemp -d "/tmp/${PACKAGE_NAME}.XXXXXX" ) || exit 1

# compile binary files
cd src
make
cd ..

# copy executable files
mkdir "${TEMPDIR}/bin"
for FILE in cookiecutter extract remove separate rm_reads; do
	cp "./src/${FILE}" "${TEMPDIR}/bin/${FILE}"
done

# copy demo and data directories
for DIR in data demo; do
	cp -R "$DIR" "${TEMPDIR}/${DIR}"
done

# adjust the path for executable files in the demo script
sed -i "s#default='.'#default='..\/bin/'#" \
	${TEMPDIR}/demo/run_demo_analysis.py

# copy LICENSE and README files
for FILE in LICENSE README.md; do
	cp "$FILE" "${TEMPDIR}/${FILE}"
done

# create a gzipped TAR archive
TEMPFILE=$( tempfile )
tar -C "/tmp" -cvf "${TEMPFILE}.tar" $( basename "$TEMPDIR" )
gzip -9 "${TEMPFILE}.tar"
mv "${TEMPFILE}.tar.gz" "$OUTPUT_FILE"

# remove compiled files
cd src
make clean
cd ..

rm -r "$TEMPDIR"

