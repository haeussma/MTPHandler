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
        +float ph*
        +Well[0..*] wells
        +float[0..*] measured_wavelengths
        +str wavelength_unit
        +AbstractSpecies[0..*] species
    }
    
    class Well {
        +InitCondition[0..*] init_conditions
        +PhotometricMeasurement[0..*] measurements
        +float volume
        +string volume_unit
        +integer x_position*
        +integer y_position*
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
        +StrictBool constant*
        +string uri
        +string creator_id
    }
    
    class AbstractSpecies {
        +string name*
        +Vessel vessel_id*
        +float init_conc
        +StrictBool constant*
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
        +Repository <sdRDM.markdown.markdownparser.MarkdownParser object at 0x130ace910>
    }
    
```