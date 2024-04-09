__version__ = '0.2.8'
import re
import os
import queue
import pandas as pd
import numpy as np
from pangu.utils import loadConfig, getLogger, BamRegionViewer, BamViewer_Error

class CoreFormater:

    def __init__( self, alignments, reference, config_yaml, impactXls, 
                  haplotypeTsv, logFile, verbose=True ):
        self.alignments = alignments
        self.reference  = reference
        self.impactXls  = impactXls
        self.hapTsv     = haplotypeTsv
        self.config     = loadConfig( config_yaml )
        self.warnings   = queue.Queue()
        self.log        = getLogger( self.__repr__(), logFile, self.warnings, stdout=verbose )
        self.region     = self.config[ 'region' ]
        self.loader     = BamRegionViewer( self.config, logger=self.log )
        self.coreVar    = self._loadCoreVar()
        self.coreMeta   = self._loadCoreMeta()

    def saveCoreFiles( self, version, outdir, backup=False ):
        '''
        Save compressed core data files to outdir.
        if backup=True, save a backup of existing files first (<file>.bak)
        '''
        for corefile, df in zip( [ 'coreVariants', 'coreMetaData' ], [ 'coreVar', 'coreMeta' ] ):
            if os.path.isfile( self.config[ corefile ] ):
                os.rename( self.config[ corefile ], self.config[ corefile ] + '.bak' )
            fname = f'{outdir}/{df.lower()}.{version}.csv.gz'
            getattr( self, df ).to_csv( fname )
            self.log.info( f'Saving core data {fname} -- UPDATE your config.yaml!' )

    def _loadCoreVar( self ):
        # load sample variants over region
        core = self.loader._makeDf( self.alignments, self.reference, self.region, qryCol='coreAllele' )
        # remove non-core alleles (alleles with dots in the name, *DD.DDD)
        core = core[~core.coreAllele.str.contains('\.')]
        # rename allele by removing the gene name
        core.coreAllele = core.coreAllele.str.extract( 'CYP2D6(\*\d+)' )
        return core
    
    def _loadCoreMeta( self ):
        import warnings
        # suppress format warning
        with warnings.catch_warnings( record=True ):
            warnings.simplefilter( "always" )
            coreDef = pd.read_excel( self.impactXls,
                                     sheet_name='Allele Function',
                                     header=1)

        alleleField   = 'Allele/cDNA/rsID'
        functionField = 'Allele Clinical Functional Status (Required)'
        
        # drop tandems
        coreDef = coreDef[ ~coreDef[ alleleField ].str.contains( 'x' ) ]
        # change core allele to integer
        coreDef[ alleleField ] = coreDef[ alleleField ].str[1:].astype( int )
        # change column name to something easier
        coreDef.rename( columns={ alleleField : 'coreAllele' }, inplace=True )
        
        coreImpacts = coreDef.set_index( 'coreAllele' )\
                             [ functionField ].map( self.config[ 'impactScores' ] )\
                             .rename( 'impact' )
        
        # make any adjustments to impact vals
        #for allele, impact in self.config[ 'impactAdjustments' ].items():
        #    coreImpacts[ allele ] = impact            
        # get count and list of var positions
        coreVars = pd.concat( [ self.coreVar.groupby( 'coreAllele' ).size().rename( 'nVars' ),
                                self.coreVar.groupby( 'coreAllele' ).pos.apply( sorted ).rename( 'varPos' ) ],
                              axis=1 )
        #change index to integers
        coreVars.index = coreVars.index.str[1:].astype(int)

        #add in rs numbers where known
        haplotypes = pd.read_csv( self.hapTsv, 
                                  sep='\t',
                                  skiprows=[0], 
                                  index_col=0 )
        core_rsID = haplotypes[ ~haplotypes.index.str.contains('\.') & haplotypes.rsID.str.contains('rs')]\
                              .dropna()\
                              .groupby( level=0 )\
                              .rsID.apply(lambda v: [*{*v}])
        core_rsID.index = core_rsID.index.str.split('\*').str[-1].astype(int)
        
        return pd.concat( [ coreImpacts, coreVars, core_rsID ], axis=1 )\
                 .fillna( dict( impact=0, nVars=0, varPos=-1 ) )\
                 .astype( dict( impact=int, nVars=int ) )

class CoreFormatter_Error(Exception):
    pass
