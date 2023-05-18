{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "b184a28a-b979-4424-bc2a-e551019558a4",
   "metadata": {},
   "outputs": [],
   "source": [
    "#!/usr/bin/env python\n",
    "\n",
    "import re\n",
    "import sys\n",
    "import json\n",
    "\n",
    "\n",
    "def tsv_parse(handle):\n",
    "    columns = handle.readline().strip('\\n').split(\"\\t\")\n",
    "    for line in handle:\n",
    "        values = line.strip('\\n').split(\"\\t\")\n",
    "        if len(columns) != len(values):\n",
    "            print(line)\n",
    "            sys.exit()\n",
    "\n",
    "\n",
    "with open(sys.arg[1]) as handle:\n",
    "    for rec in tsv_parse(handle):\n",
    "        print(json.dumps(rec))\n",
    "\n",
    "\n"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.11.3"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
