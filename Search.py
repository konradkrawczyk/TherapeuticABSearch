'''
This script 

Author: Konrad Krawczyk, NaturalAntibody

Prerequisites: ANARCI

Usage:

'''

from optparse import OptionParser
import json,os,string
from os.path import join
from anarci.anarci import run_anarci

where_i_am = os.path.dirname(os.path.abspath(__file__))

#########
#HELPERS#
#########

#based on IMGT numbering, get the CDR region
def get_region(pos,chain,scheme="imgt"):
	if scheme == 'imgt':
		if 0<=pos and pos<=26:
			return "fw"+chain+"1"
		if 27<=pos and pos<=38:
			return "cdr"+chain+"1"
		if 39<=pos and pos<=55:
			return "fw"+chain+"2"
		if 56<=pos and pos<=65:
			return "cdr"+chain+"2"
		if 66<=pos and pos<=104:
			return "fw"+chain+"3"
		if 105<=pos and pos<=117:
			return "cdr"+chain+"3"
		if 118<=pos and pos<=129:
			return "fw"+chain+"4"
	else:
		print "Scheme",scheme,'not supported'

#Number and transform the sequence into the format used in the db region->imgt->aa
#If get_germline is set to  true, germline is returned as well.
def number_and_transform(sequence,get_germline=False):

		#Get the numbered sequence out...
		numbered_sequence = number_sequence(sequence)
		#Check that the input was ok by virtue of being able to number it ok
		if numbered_sequence[1][0] == None:
			print "There was a problem numbering sequence",sequence,'please review your input'
			quit()
		germlines = {'V':'','J':'','species':''}
		if get_germline == True:
			germlines['V'] = numbered_sequence[2][0][0]['germlines']['v_gene'][0][1]
			germlines['J'] = numbered_sequence[2][0][0]['germlines']['j_gene'][0][1]
			germlines['species'] = numbered_sequence[2][0][0]['germlines']['v_gene'][0][0]
		if numbered_sequence[1][0] == None:
			return None
		#Extract the chain.
		chain = numbered_sequence[2][0][0]['chain_type'].lower()
		if chain == 'k':
			chain = 'l'
		#Make it sortable and identifiable.
		numbered_sequence = transform_anarci_output(numbered_sequence[1][0][0][0])
		
		#transform to output as seen in the db
		organized = {}
		for elem in numbered_sequence:
			#We do not care for gaps.
			if numbered_sequence[elem] == '-':
				continue
			region = get_region(int(elem[0]),chain)
			if region not in organized:
				organized[region] = {}
			 
			organized[region][str(elem[0])+elem[1].replace(' ','')] = numbered_sequence[elem]
		if get_germline == True:
			return {'numbered':organized,'germlines':germlines}
		else:
			return organized

#Number sequence using anarci
def number_sequence(query_seq,force_human = False,assign_germline=True,scheme="imgt"):
	#Number the query sequence
	#
	if force_human == True:
		res = run_anarci([('q',query_seq)],scheme=scheme,assign_germline=True,allowed_species=["human"])
	else:
		res = run_anarci([('q',query_seq)],scheme=scheme,assign_germline=True)
	
	return res

#Transform anarci list into dictionary.
def transform_anarci_output(anarci_output):
		
	sequence = {}
	for elem in anarci_output:
		
		num,insertion = elem[0]
		insertion = insertion.replace(' ','')
		if elem[1] == '-':#We do not care for gaps.
			continue
		sequence[(num,insertion)] = elem[1]
		
	return sequence

#Sequence identity on the region-annotated format.
#TODO introduce efficient representation search.
#If normalization is turned on, identity will be averaged.
def sequence_identity(query,template,normalize=False):
		
	total = 0
	matches = 0
	regions = {}
	
	for region in query:
		#Calculate alignment to regions.
		r_total = 0
		regions[region] = {'s':{},'id':0}
		for imgt in query[region]:
			#print imgt,total
			r_total+=1
			total+=1
			try:
				regions[region]['s'][imgt] = template[region][imgt]
				if query[region][imgt] == template[region][imgt]:
					matches+=1
					regions[region]['id']+=1
			
			except KeyError:#Region or imgt entry are not in the template.
				pass
		#Region sequence identity
		regions[region]['id'] = normalize_number(100*(float(regions[region]['id'])/float(r_total)))
	
	sid = normalize_number(100*(float(matches)/float(total)))	

	#Are we missing any regions?
	missing_regions = []
	if normalize == True:

		sid_2 = sequence_identity(template,query,normalize=False)
		
		for region in regions:
			try:
				sid_region = float(regions[region]['id']+sid_2['regions'][region]['id'])/2.0
				
				regions[region]['id'] = normalize_number(sid_region)
			except KeyError:
				missing_regions.append(region)
		sid = float(sid+sid_2['sid'])/2.0
		
	
	return {'sid':normalize_number(sid),'regions':regions,'missing':missing_regions}

#Drop the number to two floating points.
def normalize_number(num):
	return int((100*num))/100.0

#Align two sequences and print the results onto the command line.
def align_sequences(s1,s2,show=False):
	
	
	#Calculate the sequence identity - average between t-q and t-q	
	identity_1 = sequence_identity(s1,s2)['sid']
	identity_2 = sequence_identity(s2,s1)['sid']

	identity = int((identity_1+identity_2)/2.0)

	s1 = get_sorted_sequence(s1)
	s2 = get_sorted_sequence(s2)
	
	#Get the keys.
	ks = []
	for imgt in s1:
		ks.append(imgt)
	for imgt in s2:
		if imgt not in ks:
			ks.append(imgt)
	s1_out = ''
	s2_out = ''
	ali = ''
	region_data = ''
	for imgt in sorted(ks):
		#The alignment quality.
		if imgt in s1 and imgt in s2:
			if s1[imgt] == '-' and s2[imgt] == '-':
				continue
			if s1[imgt] == s2[imgt]:
				ali+='|'
			else:
				ali+='.'
		else:
			ali+='.'
		#Amino acid identities.
		if imgt in s1:
			s1_out+=s1[imgt]
		else:
			s1_out+='-'
		if imgt in s2:
			s2_out+=s2[imgt]
		else:
			s2_out+='-'
		#Region annotation#TODO hardcoded this bit.
		reg = get_region(imgt[0],'H')
		if 'cdr' in reg:
			region_data+='^'
		else:
			region_data+=' '

	if show ==True:
		print s1_out
		print ali		
		print s2_out
		print region_data
		print "Identity=",identity,'%\n\n'
	else:
		out = s1_out+'\n'+ali+'\n'+s2_out+'\n'+region_data
		return out

#Transrofm given json-formatted sequence to proper sequence file.
def get_sorted_sequence(data,target_region=None):

	#We use imgt numbering. this means that certain entries go after insertions
	#32 31A 31B ... 33A 33B 33
	#60 60A 60B ... 61A 61B 61
	#111 111A 111B ... 112A 112B 112C
	unusual_inserts = [33,61,112]
	
	sequence = read_sequence(data,target_region=target_region)
	
	p_sequence = ""
	#Dealing with IMGT ordering.
	for uid in unusual_inserts:
		
		imgt_id = (uid,'','')
		
		if imgt_id in sequence:
			aa = sequence[imgt_id]
			del sequence[imgt_id]
			new_id = (uid,1000,'')
			sequence[new_id] = aa	
		for alpha in string.ascii_uppercase:
			imgt_id = (uid,'',alpha)
			if imgt_id in sequence:
				aa = sequence[imgt_id]
				del sequence[imgt_id]
				new_id = (uid,1000-ord(alpha),alpha)
				sequence[new_id] = aa
	
	return sequence

#Read the sequence format and return python-formatted variety
def read_sequence(data,target_region=None):
	sequence = {}
	for region in data:
		if target_region != None and region!=target_region:
			continue
		for imgt_id in data[region]:
			
			#See if we have insertions?
			try:
				_id = int(imgt_id)
				_id = (_id,'','')
			except ValueError:
				insertion = imgt_id[-1]
				_id = int(imgt_id[0:(len(imgt_id)-1)])
				_id = (_id,'',insertion)
			sequence[_id] = data[region][imgt_id]

	return sequence

#This is the main class that is used for searching - for multiple searches it is advisable to re-use the object as it caches the therapetuics database in memory
class Search:

	def __init__(self):
		
		#Pre-load the set of therapetuics.
		self.therapeutics = json.load(open(join('data','numbered.json')))
		#print self.therapeutics

	#Perform the search given parameters formatted as following:
	#parameters = { 'heavy_sequence':"EVQL...",#sequence of the heavy chain
	#		'light_sequence':"DIV...",#sequence of the light chain
	#		'cutoff':80,#sequence identity cutoff
	#		}
	
	#Search for a single sequence - needs to be numbered
	def single_search(self,numbered=None,cutoff=80):
		
		#Go through therapeutics and perform search.
		#These are the sequence matches
		sequence_matches = []
		for therapeutic in self.therapeutics:
			for chain_id in ['1','2']:#Account for bispecifics
				if self.therapeutics[therapeutic][numbered['c']+chain_id+'_numbered']==None:
					continue
				sid_results = sequence_identity(numbered['n'],self.therapeutics[therapeutic][numbered['c']+chain_id+'_numbered'])
				#Check if the antibody sequence checks the cutoff box.
				SEQID = None
				if parameters['region'] == 'V':
					SEQID = sid_results['sid']#Default - v region identity
				else:
					for region in sid_results['regions']:
						if 'fw' in region:
							r = 'fw'+region[-1]
						
						elif 'cdr' in region:
						
							r = 'cdr'+region[-1]
						if r == parameters['region']:
							SEQID = float(sid_results['regions'][region]['id'])
							
				if SEQID == None:
					continue #Means that there is no region like that.
				if SEQID >= cutoff:
					match = {'name':therapeutic+'-'+numbered['c']+chain_id,'regions':sid_results['regions'],'target':self.therapeutics[therapeutic]['target'],'numbered':self.therapeutics[therapeutic][numbered['c']+chain_id+'_numbered']}
					sequence_matches .append( (SEQID,match) )

		return sequence_matches

	#Perform the search.
	def perform_search(self,parameters):
		
		#These are the formatted results
		results = {
				'heavy_sequence_results':None,#Results for heavy chain only
				'light_sequence_results':None,#Results for light chain only
				'combined_results':None#Results for heavy and light chain combined
			}

		#Number the sequence.
		for s in ['heavy_sequence','light_sequence']:
			if parameters[s]==None:#Are we performing search on just one sequence?
				
				continue
			numbered = number_and_transform(parameters[s],get_germline=True)
			chain_type = numbered['germlines']['V'][2]
			if chain_type == "K":
				chain_type = "L"
			numbered = {'n':numbered['numbered'],'c':chain_type}
			parameters[s+'_numbered'] = numbered
			results[s+'_results'] = self.single_search(numbered=parameters[s+'_numbered'],cutoff=parameters['cutoff'])
			

		#Digest the results
		self.digest_results(results,parameters)

	#Digest the results.
	def digest_results(self,results,parameters):
		for s in ['heavy_sequence','light_sequence']:
			if results[s+'_results']!=None:

				self.write_out_results(results[s+'_results'],join(parameters['output'],s),parameters[s+'_numbered'])
		
		#Process combined results - this is for the V region only
		if results['heavy_sequence_results']!=None and results['light_sequence_results']!=None and parameters['region'] == 'V':
			#Collect the antibodies which tick the cutoff box on both counts.
			heavy_names = {}
			light_names = {}
			for ab in sorted(results['heavy_sequence_results']):
				sid = ab[0]
				name = ab[1]['name'].split('-')[0]
				heavy_names[name] = sid
			for ab in sorted(results['light_sequence_results']):
				sid = ab[0]
				name = ab[1]['name'].split('-')[0]
				light_names[name] = sid
			f = open(join(parameters['output'],'combined.tsv'),'w')
			f.write('Therapeutic name\tHeavy chain V Identity (%)\tLight chain V Identity (%)')
			for name in heavy_names:
				if name in light_names:
					f.write(name+'\t'+str(heavy_names[name])+'\t'+str(light_names[name])+'\n')

			f.close()
		
				
	#shorthand for writing out a single piece of results
	def write_out_results(self,results,out_name,query):
		out_identities = open(out_name+'_identities.tsv','w')
		out_alignments = open(out_name+'_alignments.txt','w')
		out_identities.write('Therapeutic Name\tTarget\tV identity(%)\tCDR1 identity (%)\tCDR2 identity (%)\tCDR2 identity (%)\tFW1 identity (%)\tFW2 identity (%)\tFW3 identity (%)\tFW4 identity (%)\n')
		for elem in sorted(results):
			
			v_sequence_identity = str(normalize_number(elem[0]))
			data = elem[1]
			therapeutic_name = data['name']
			target = data['target'].decode('utf-8')
			#Alignments
			out_alignments.write('Therapeutic = '+therapeutic_name+'\n')
			out_alignments.write(align_sequences(query['n'],data['numbered'])+'\n')
			
			
			#Region identities
			region_ids= ""
			for region in sorted(data['regions']):
				region_ids+='\t'+str(normalize_number(data['regions'][region]['id']))
			out_identities.write(therapeutic_name+'\t'+target+'\t'+v_sequence_identity+region_ids+'\n')
		out_alignments.close()
		out_identities.close()

if __name__ == '__main__':
	

	usage = '''
python Search.py  --heavy-sequence QVQLQQSGSELKKPGASVKVSCKASGYTFTNYGMNWVKQAPGQGLKWMGWINTYTGEPTYTDDFKGRFAFSLDTSVSTAYLQISSLKADDTAVYFCARGGFGSSYWYFDVWGQGSLVTVSS --light-sequence DIQMTQSPSSLSASVGDRVTITCSASQDISNYLNWYQQKPGKAPKVLIYFTSSLHSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYSTVPWTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCL	
	'''

	parser = OptionParser()
	parser.add_option("--heavy-sequence", dest="heavy_sequence",
                  help="Heavy Chain Sequence",)
	parser.add_option("--light-sequence", dest="light_sequence",
                  help="Light Chain Sequence")
	parser.add_option("--cutoff", dest="cutoff",default=80,
                  help="Sequence identity cutoff")
	parser.add_option("--region", dest="region",default="V",
                  help="Region to perform the search on V [default], cdr1,cdr2,cdr3,fw1,fw2,fw3,fw4")
	parser.add_option("--output", dest="output",default="search_out",
                  help="Directory to write the output files into")
	
	(options, args) = parser.parse_args()

	#Check input, we need at least one heavy and one light chain.
	if not options.heavy_sequence and not options.light_sequence:
		print "Your input is missing the sequence information, use this code as follows:"
		print usage
		quit()

	
	#Collect input.
	parameters = {
			'heavy_sequence':None,
			'light_sequence':None,
			'cutoff':None,
			'output':'search_out',
			}
	#Output directory
	if not os.path.exists(options.output):

		os.mkdir(options.output)
	#Heavy chain
	if options.heavy_sequence:
		parameters["heavy_sequence"] = options.heavy_sequence
	#Light chain
	if options.light_sequence:
		parameters["light_sequence"] = options.light_sequence
	#If not user-specified will go to the default location.
	parameters['output'] = options.output
	#If not user-specified will go to the default parameter.
	parameters['cutoff'] = float(options.cutoff)
	#If not user-specified will go to the default parameter.
	parameters['region'] = options.region
	allowed_regions = ['V','cdr1','cdr2','cdr3','fw1','fw2','fw3','fw4']
	if parameters['region'] not in allowed_regions:
		print "You specified region that is disallowed:",options.region
		print "Allowed regions to perform search on are: ",allowed_regions
	print "Results will be written to",options.output
	

	#Define the search object
	search = Search()
	
	search.perform_search(parameters)


	
