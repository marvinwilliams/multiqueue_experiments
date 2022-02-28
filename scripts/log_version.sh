#!/bin/bash

git -C .. diff HEAD > experiments.patch
git -C ../multiqueue diff HEAD > multiqueue.patch
