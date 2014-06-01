#!/usr/bin/python

# 08.05.2014 Lluis Revilla
# This script parses hmmer3 text files from the website 
# (without using biopython 1.61 or above and its parser.)
# Given a hmmer file from hmmer.janelia.org select the better e-value for each alignment and the non overlapping domains for each query
# Requires the modules: re, argparse, and csv
# It is thought to run from command line

import re
import argparse # Manage arguments
import csv # To create the csv file. 
parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter, description = "This script given a text file from hmmscan (without any formatting option) select the best alignment (best e-value) for each domain and then outputs the non overlapping domains for each protein in a csv file." )
parser.add_argument("hmmer_file", help = """The output of the hmmscan runned locally over a multifasta protein file.""", type=str)
parser.add_argument('output_file', help = """Name of the csv file comma separated with the selected domains.""", type = str )
parser.add_argument('-threshold', help = """If you want to include domains under the threshold, use this option. The name of these domains are marked with a "!" at the end.""", action = 'store_false')
parser.add_argument('-version', action = 'version', version = '1')
args = parser.parse_args()

########################
## Open and read file ##
########################

f = open(args.hmmer_file, 'r')
document = f.readlines()
f.close()

####################################
## Search and store report values ##
####################################
# The regular expressions to select the lines
query_regex = re.compile("^(Query:)")
domain_regex = re.compile("^ {3,10}((-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*){8}.+")
inclusion_regex = re.compile("^(  ------ inclusion threshold ------\n)")
alignments_regex = re.compile("[?!]")

data = {}
names_query = []
for line in document:
	# Skip the comments lines to make it faster directly continue
	if line.startswith("#") or '*' in line:
		continue
	elif re.findall(query_regex, line) != []:
		threshold = False
		query = line.split()[1] # Protein name
		names_query.append(query) # To get an ordered set of names
		data[query] = {}

		names_domain = []
		n_domain = 0
		n_hits = 0
		n_domains = 0
		data[query]["n_hits"] = n_hits
	elif re.findall(domain_regex, line) != []:
		domain = line.split()[8] # Name of domain
		if threshold == True:
			if args.threshold == False:
				domain += "!"
			else:
				continue # Skip the domains under threshold
		names_domain.append(domain)
		data[query][domain] = {}
		n_hits += 1
		data[query]["n_hits"] +=1 
		description = " ".join(line.split()[9:]) # Description domain
		data[query][domain]["Description"] = description
		data[query][domain]["e_value"] = float(line.split()[0])
	elif re.findall(inclusion_regex, line) != []:
		threshold = True
	elif re.findall(alignments_regex, line) != []:
		lines = line.split()
		if n_domains > n_hits:
			continue
		if int(lines[0]) == 1: #It is the first alignment
			n_domains += 1
			if n_domains > n_hits:
				continue
			if names_domain != [] and n_domain < len(names_domain):
				domain = names_domain[n_domain]
			else:
				continue #No more domains 
			try:
				data[query][domain]["e_value_alignment"] > float(lines[5])
			except:
				alignment = data[query][domain]
				alignment["e_value_alignment"] = float(lines[5])
				alignment["hit_start"] = int(lines[6])
				alignment["hit_end"]  = int(lines[7])
				alignment["query_start"]  = query_start = int(lines[9])
				alignment["query_end"]  = query_end = int(lines[10])
				alignment["Length"]  = query_end - query_start
			if n_domain < len(names_domain):
				domain = names_domain[n_domain]
				n_domain += 1
		elif alignment["e_value_alignment"] > float(lines[5]):
			alignment["e_value_alignment"] = float(lines[5])
			alignment["hit_start"] = int(lines[6])
			alignment["hit_end"]  = int(lines[7])
			alignment["query_start"]  = query_start = int(lines[9])
			alignment["query_end"]  = query_end = int(lines[10])
			alignment["Length"]  = query_end - query_start

###########################################
## Select the domains that don't overlap ##
###########################################
selected = {}
for query in names_query:
	lengths = [] # Used to select the values
	domains_all = [] #Store the list of hits

	# Remove the n_hits from the name of domains...
	for domain in data[query]:
		if domain == "n_hits" or domain == "e_value":
			pass
		else:
			domains_all.append(domain)
		
	if domains_all == []: #Skip empty queries
		continue
	domains = [] # Store the desired values
	from_positions = [] # Used to select the values
	to_positions = [] # Used to select the values

	for domain in domains_all:
		try:
			lengths.append(data[query][domain]["Length"])
		except:
			print data[query]

	lengths = sorted(lengths)[::-1] #Set the longest in the first position [0]

	for length in lengths:
		for domain in domains_all:
			info = data[query][domain]
			# Set values to compare between the hits
			from_position = info["query_start"]
			to_position = info["query_end"]

			if length == info["Length"]:
				# The first to add is the longest one
				if len(domains) == 0:
					from_positions.append(from_position)
					to_positions.append(to_position)
					domains.append(domain)
					break
				# If there is one already added
				elif len(domains) >= 1:
					min_from = min(from_positions)
					max_to = max(to_positions)
					# Outside the domains already added
					if to_position < min_from or from_position > max_to:
						from_positions.append(from_position)
						to_positions.append(to_position)
						domains.append(domain)
						break 
					# The domain could be between two domains
					else:
						# It is not added in from_position sorted
						# Sorting it we ensure to 
						from_positions_1 = sorted(from_positions)
						to_positions_1 = sorted(to_positions)
						for i in range(len(from_positions_1)-1):
							if from_position > to_positions_1[i] and to_position < from_positions_1[1+i]:
								from_positions.append(from_position)
								to_positions.append(to_position)
								domains.append(domain)
								break
							else:
								pass #Didn't fit, it is overlapping a % of it. A parameter to accept some of them?

	# Store the selected domains
	selected[query] = domains 

#######################
## Print out the csv ##
#######################
with open(args.output_file, "w+") as csvfile:
	a = csv.writer(csvfile, delimiter=',')
	first_row = ["Query name", "Description", "Name",
	 "significance","hits", "hit_start", "hit_end", "query_start",
	 "query_end", "Length"]
	a.writerow(first_row)

	for query in names_query:
		if query not in selected.keys():
			continue #Skip the empty ones
		for domain in selected[query]:
			info = data[query][domain]
			data_row = [query,
			 info["Description"],
			 domain,
			 data[query][domain]["e_value"],
			 data[query]["n_hits"], # Number of hits
			 info["hit_start"],
			 info["hit_end"],
			 info["query_start"],
			 info["query_end"],
			 info["Length"]]
			a.writerow(data_row)
