```mermaid
classDiagram
    AbstractSpecies <-- Protein
    AbstractSpecies <-- Complex
    AbstractSpecies <-- Reactant
    EnzymeMLDocument *-- Creator
    EnzymeMLDocument *-- Vessel
    EnzymeMLDocument *-- Protein
    EnzymeMLDocument *-- Complex
    EnzymeMLDocument *-- Reactant
    EnzymeMLDocument *-- Reaction
    EnzymeMLDocument *-- KineticParameter
    EnzymeMLDocument *-- Measurement
    EnzymeMLDocument *-- File
    AbstractSpecies *-- Vessel
    Protein *-- SBOTerm
    Complex *-- SBOTerm
    Reactant *-- SBOTerm
    Reaction *-- SBOTerm
    Reaction *-- ReactionElement
    Reaction *-- KineticModel
    ReactionElement *-- SBOTerm
    ReactionElement *-- AbstractSpecies
    KineticModel *-- SBOTerm
    KineticModel *-- KineticParameter
    KineticParameter *-- SBOTerm
    Measurement *-- MeasurementData
    MeasurementData *-- AbstractSpecies
    MeasurementData *-- Replicate
    Replicate *-- DataTypes
    Replicate *-- AbstractSpecies
    Plate *-- Well
    Plate *-- AbstractSpecies
    Well *-- PhotometricMeasurement
    Well *-- InitCondition
    PhotometricMeasurement *-- BlankState
    InitCondition *-- AbstractSpecies
    BlankState *-- AbstractSpecies
    AbstractSpecies *-- Vessel
    
    class Plate {
        +integer n_rows*
        +integer n_columns*
        +datetime date_measured
        +float[0..*] times
        +str time_unit
        +float[0..*] temperatures*
        +str temperature_unit*
        +float max_volume
        +str max_volume_unit
        +Well[0..*] wells
        +float[0..*] measured_wavelengths
        +str wavelength_unit
        +AbstractSpecies[0..*] species
    }
    
    class Well {
        +float ph*
        +integer x_position*
        +integer y_position*
        +InitCondition[0..*] init_conditions
        +PhotometricMeasurement[0..*] measurements
        +float volume
        +string volume_unit
    }
    
    class PhotometricMeasurement {
        +float wavelength*
        +str wavelength_unit*
        +float[0..*] absorptions*
        +BlankState[0..*] blank_states
    }
    
    class InitCondition {
        +AbstractSpecies species_id*
        +float init_conc*
        +str conc_unit*
    }
    
    class BlankState {
        +AbstractSpecies species_id*
        +bool contributes_to_signal*
    }
    
    class Vessel {
        +string name*
        +posfloat volume*
        +string unit*
        +boolean constant*
        +string uri
        +string creator_id
    }
    
    class AbstractSpecies {
        +string name*
        +Vessel vessel_id*
        +float init_conc
        +boolean constant*
        +string unit
        +string uri
        +string creator_id
    }
    
    class SBOTerm {
        << Enumeration >>
        +BIOCHEMICAL_REACTION
        +ACID_BASE_REACTION
        +CONFORMATIONAL_TRANSITION
        +CONVERSION
        +DEGRADATION
        +DISSOCIATION
        +IONISATION
        +ISOMERISATION
        +NON_COVALENT_BINDING
        +REDOX_REACTION
        +SPONTANEOUS_REACTION
        +PROTEIN
        +GENE
        +SMALL_MOLECULE
        +ION
        +RADICAL
        +INTERACTOR
        +SUBSTRATE
        +PRODUCT
        +CATALYST
        +INHIBITOR
        +ESSENTIAL_ACTIVATOR
        +NON_ESSENTIAL_ACTIVATOR
        +POTENTIATOR
        +MACROMOLECULAR_COMPLEX
        +PROTEIN_COMPLEX
        +DIMER
        +MICHAELIS_MENTEN
        +K_CAT
        +K_M
        +V_MAX
    }
    
    class DataTypes {
        << Enumeration >>
        +CONCENTRATION
        +ABSORPTION
        +FEED
        +BIOMASS
        +CONVERSION
        +PEAK_AREA
    }
    
    class https://github.com/EnzymeML/enzymeml-specifications.git {
        << External Object >>
        +Repository objects=[{'name': 'Vessel', 'docstring': 'This object describes vessels in which the experiment has been carried out. These can include any type of vessel used in biocatalytic experiments.', 'attributes': [{'name': 'name', 'required': True, 'type': ['string'], 'description': 'Name of the used vessel.', 'template_alias': 'Name'}, {'name': 'volume', 'required': True, 'type': ['posfloat'], 'description': 'Volumetric value of the vessel.', 'template_alias': 'Volume value'}, {'name': 'unit', 'required': True, 'type': ['string'], 'description': 'Volumetric unit of the vessel.', 'template_alias': 'Volume unit'}, {'name': 'constant', 'required': True, 'type': ['boolean'], 'description': 'Whether the volume of the vessel is constant or not.', 'default': 'True'}, {'name': 'uri', 'required': False, 'default': None, 'type': ['string'], 'description': 'URI of the vessel.'}, {'name': 'creator_id', 'required': False, 'default': None, 'type': ['string'], 'description': 'Unique identifier of the author.'}], 'type': 'object', 'subtypes': [], 'module': 'species'}, {'name': 'AbstractSpecies', 'docstring': 'This object is used to inherit basic attributes common to all species used in the data model.', 'attributes': [{'name': 'name', 'required': True, 'type': ['string'], 'description': 'None'}, {'name': 'vessel_id', 'required': True, 'reference': 'Vessel.id', 'type': ['Vessel'], 'description': 'None'}, {'name': 'init_conc', 'required': False, 'default': None, 'type': ['float'], 'description': 'None'}, {'name': 'constant', 'required': True, 'type': ['boolean'], 'description': 'None'}, {'name': 'unit', 'required': False, 'default': None, 'type': ['string'], 'description': 'None'}, {'name': 'uri', 'required': False, 'default': None, 'type': ['string'], 'description': 'None'}, {'name': 'creator_id', 'required': False, 'default': None, 'type': ['string'], 'description': 'None'}], 'type': 'object', 'subtypes': [], 'module': 'species'}] enums=[{'name': 'SBOTerm', 'docstring': 'These are a small fraction of the SBOTerms defined for the SBML markup language.', 'mappings': [{'key': 'BIOCHEMICAL_REACTION', 'value': '"SBO:0000176"'}, {'key': 'ACID_BASE_REACTION', 'value': '"SBO:0000208"'}, {'key': 'CONFORMATIONAL_TRANSITION', 'value': '"SBO:0000181"'}, {'key': 'CONVERSION', 'value': '"SBO:0000182"'}, {'key': 'DEGRADATION', 'value': '"SBO:0000179"'}, {'key': 'DISSOCIATION', 'value': '"SBO:0000180"'}, {'key': 'IONISATION', 'value': '"SBO:0000209"'}, {'key': 'ISOMERISATION', 'value': '"SBO:0000377"'}, {'key': 'NON_COVALENT_BINDING', 'value': '"SBO:0000177"'}, {'key': 'REDOX_REACTION', 'value': '"SBO:0000200"'}, {'key': 'SPONTANEOUS_REACTION', 'value': '"SBO:0000672"'}, {'key': 'PROTEIN', 'value': '"SBO:0000252"'}, {'key': 'GENE', 'value': '"SBO:0000251"'}, {'key': 'SMALL_MOLECULE', 'value': '"SBO:0000247"'}, {'key': 'ION', 'value': '"SBO:0000327"'}, {'key': 'RADICAL', 'value': '"SBO:0000328"'}, {'key': 'INTERACTOR', 'value': '"SBO:0000336"'}, {'key': 'SUBSTRATE', 'value': '"SBO:0000015"'}, {'key': 'PRODUCT', 'value': '"SBO:0000011"'}, {'key': 'CATALYST', 'value': '"SBO:0000013"'}, {'key': 'INHIBITOR', 'value': '"SBO:0000020"'}, {'key': 'ESSENTIAL_ACTIVATOR', 'value': '"SBO:0000461"'}, {'key': 'NON_ESSENTIAL_ACTIVATOR', 'value': '"SBO:0000462"'}, {'key': 'POTENTIATOR', 'value': '"SBO:0000021"'}, {'key': 'MACROMOLECULAR_COMPLEX', 'value': '"SBO:0000296"'}, {'key': 'PROTEIN_COMPLEX', 'value': '"SBO:0000297"'}, {'key': 'DIMER', 'value': '"SBO:0000607"'}, {'key': 'MICHAELIS_MENTEN', 'value': '"SBO:0000028"'}, {'key': 'K_CAT', 'value': '"SBO:0000025"'}, {'key': 'K_M', 'value': '"SBO:0000027"'}, {'key': 'V_MAX', 'value': '"SBO:0000186"'}], 'type': 'enum'}, {'name': 'DataTypes', 'docstring': 'These values are used to determine the type of time course data.', 'mappings': [{'key': 'CONCENTRATION', 'value': '"conc"'}, {'key': 'ABSORPTION', 'value': '"abs"'}, {'key': 'FEED', 'value': '"feed"'}, {'key': 'BIOMASS', 'value': '"biomass"'}, {'key': 'CONVERSION', 'value': '"conversion"'}, {'key': 'PEAK_AREA', 'value': '"peak-area"'}], 'type': 'enum'}] inherits=[{'parent': 'AbstractSpecies', 'child': 'Protein'}, {'parent': 'AbstractSpecies', 'child': 'Complex'}, {'parent': 'AbstractSpecies', 'child': 'Reactant'}] compositions=[{'container': 'EnzymeMLDocument', 'module': 'Creator'}, {'container': 'EnzymeMLDocument', 'module': 'Vessel'}, {'container': 'EnzymeMLDocument', 'module': 'Protein'}, {'container': 'EnzymeMLDocument', 'module': 'Complex'}, {'container': 'EnzymeMLDocument', 'module': 'Reactant'}, {'container': 'EnzymeMLDocument', 'module': 'Reaction'}, {'container': 'EnzymeMLDocument', 'module': 'KineticParameter'}, {'container': 'EnzymeMLDocument', 'module': 'Measurement'}, {'container': 'EnzymeMLDocument', 'module': 'File'}, {'container': 'AbstractSpecies', 'module': 'Vessel'}, {'container': 'Protein', 'module': 'SBOTerm'}, {'container': 'Complex', 'module': 'SBOTerm'}, {'container': 'Reactant', 'module': 'SBOTerm'}, {'container': 'Reaction', 'module': 'SBOTerm'}, {'container': 'Reaction', 'module': 'ReactionElement'}, {'container': 'Reaction', 'module': 'KineticModel'}, {'container': 'ReactionElement', 'module': 'SBOTerm'}, {'container': 'ReactionElement', 'module': 'AbstractSpecies'}, {'container': 'KineticModel', 'module': 'SBOTerm'}, {'container': 'KineticModel', 'module': 'KineticParameter'}, {'container': 'KineticParameter', 'module': 'SBOTerm'}, {'container': 'Measurement', 'module': 'MeasurementData'}, {'container': 'MeasurementData', 'module': 'AbstractSpecies'}, {'container': 'MeasurementData', 'module': 'Replicate'}, {'container': 'Replicate', 'module': 'DataTypes'}, {'container': 'Replicate', 'module': 'AbstractSpecies'}] external_objects={} namespaces={} add_id_field=True
    }
    
```