#!/bin/bash

cat map.commands.txt | xargs -I cathy -P10 sh -c cathy
