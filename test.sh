#!/bin/bash

set -e

source envvar.sh

source py27env/bin/activate

cd slingen
python temp_main.py

deactivate
