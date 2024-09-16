## This is a generated file. Do not modify it manually!

from __future__ import annotations

from enum import Enum
from typing import Generic, Optional, TypeVar
from uuid import uuid4

from pydantic import BaseModel, ConfigDict, Field

# Filter Wrapper definition used to filter a list of objects
# based on their attributes
Cls = TypeVar("Cls")


class FilterWrapper(Generic[Cls]):
    """Wrapper class to filter a list of objects based on their attributes"""

    def __init__(self, collection: list[Cls], **kwargs):
        self.collection = collection
        self.kwargs = kwargs

    def filter(self) -> list[Cls]:
        for key, value in self.kwargs.items():
            self.collection = [
                item for item in self.collection if self._fetch_attr(key, item) == value
            ]
        return self.collection

    def _fetch_attr(self, name: str, item: Cls):
        try:
            return getattr(item, name)
        except AttributeError:
            raise AttributeError(f"{item} does not have attribute {name}")


# JSON-LD Helper Functions
def add_namespace(obj, prefix: str | None, iri: str | None):
    """Adds a namespace to the JSON-LD context

    Args:
        prefix (str): The prefix to add
        iri (str): The IRI to add
    """
    if prefix is None and iri is None:
        return
    elif prefix and iri is None:
        raise ValueError("If prefix is provided, iri must also be provided")
    elif iri and prefix is None:
        raise ValueError("If iri is provided, prefix must also be provided")

    obj.ld_context[prefix] = iri  # type: ignore


def validate_prefix(term: str | dict, prefix: str):
    """Validates that a term is prefixed with a given prefix

    Args:
        term (str): The term to validate
        prefix (str): The prefix to validate against

    Returns:
        bool: True if the term is prefixed with the prefix, False otherwise
    """

    if isinstance(term, dict) and not term["@id"].startswith(prefix + ":"):
        raise ValueError(f"Term {term} is not prefixed with {prefix}")
    elif isinstance(term, str) and not term.startswith(prefix + ":"):
        raise ValueError(f"Term {term} is not prefixed with {prefix}")


# Model Definitions


class Plate(BaseModel):
    model_config: ConfigDict = ConfigDict(  # type: ignore
        validate_assigment=True,
    )  # type: ignore

    temperatures: list[float]
    temperature_unit: UnitDefinition
    id: Optional[str] = Field(default=None)
    name: Optional[str] = Field(default=None)
    wells: list[Well] = Field(default_factory=list)
    species: list[Species | Molecule | Protein] = Field(default_factory=list)
    date_measured: Optional[str] = Field(default=None)
    times: list[float] = Field(default_factory=list)
    time_unit: Optional[UnitDefinition] = Field(default=None)

    # JSON-LD fields
    ld_id: str = Field(
        serialization_alias="@id", default_factory=lambda: "md:Plate/" + str(uuid4())
    )
    ld_type: list[str] = Field(
        serialization_alias="@type",
        default_factory=lambda: [
            "md:Plate",
        ],
    )
    ld_context: dict[str, str | dict] = Field(
        serialization_alias="@context",
        default_factory=lambda: {
            "md": "https://github.com/FAIRChemistry/MTPHandler",
        },
    )

    def filter_wells(self, **kwargs) -> list[Well]:
        """Filters the wells attribute based on the given kwargs

        Args:
            **kwargs: The attributes to filter by.

        Returns:
            list[Well]: The filtered list of Well objects
        """

        return FilterWrapper[Well](self.wells, **kwargs).filter()

    def filter_species(self, **kwargs) -> list[Species]:
        """Filters the species attribute based on the given kwargs

        Args:
            **kwargs: The attributes to filter by.

        Returns:
            list[Species]: The filtered list of Species objects
        """

        return FilterWrapper[Species](self.species, **kwargs).filter()

    def set_attr_term(
        self,
        attr: str,
        term: str | dict,
        prefix: str | None = None,
        iri: str | None = None,
    ):
        """Sets the term for a given attribute in the JSON-LD object

        Example:
            # Using an IRI term
            >> obj.set_attr_term("name", "http://schema.org/givenName")

            # Using a prefix and term
            >> obj.set_attr_term("name", "schema:givenName", "schema", "http://schema.org")

            # Usinng a dictionary term
            >> obj.set_attr_term("name", {"@id": "http://schema.org/givenName", "@type": "@id"})

        Args:
            attr (str): The attribute to set the term for
            term (str | dict): The term to set for the attribute

        Raises:
            AssertionError: If the attribute is not found in the model
        """

        assert (
            attr in self.model_fields
        ), f"Attribute {attr} not found in {self.__class__.__name__}"

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_context[attr] = term

    def add_type_term(
        self, term: str, prefix: str | None = None, iri: str | None = None
    ):
        """Adds a term to the @type field of the JSON-LD object

        Example:
            # Using a term
            >> obj.add_type_term("https://schema.org/Person")

            # Using a prefixed term
            >> obj.add_type_term("schema:Person", "schema", "https://schema.org/Person")

        Args:
            term (str): The term to add to the @type field
            prefix (str, optional): The prefix to use for the term. Defaults to None.
            iri (str, optional): The IRI to use for the term prefix. Defaults to None.

        Raises:
            ValueError: If prefix is provided but iri is not
            ValueError: If iri is provided but prefix is not
        """

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_type.append(term)

    def add_to_wells(
        self,
        id: str,
        x_pos: int,
        y_pos: int,
        ph: Optional[float] = None,
        init_conditions: list[InitCondition] = [],
        measurements: list[PhotometricMeasurement] = [],
        volume: Optional[float] = None,
        volume_unit: Optional[UnitDefinition] = None,
        **kwargs,
    ):
        params = {
            "id": id,
            "x_pos": x_pos,
            "y_pos": y_pos,
            "ph": ph,
            "init_conditions": init_conditions,
            "measurements": measurements,
            "volume": volume,
            "volume_unit": volume_unit,
        }

        if "id" in kwargs:
            params["id"] = kwargs["id"]

        self.wells.append(Well(**params))

        return self.wells[-1]

    def add_species_to_species(
        self,
        ld_id: Optional[str] = None,
        id: Optional[str] = None,
        name: Optional[str] = None,
    ):
        params = {
            "ld_id": ld_id,
            "id": id,
            "name": name,
        }
        if not ld_id:
            params.pop("ld_id")

        species = Species(**params)

        self.species.append(species)

        return species

    def add_molecule_to_species(
        self,
        smiles: Optional[str] = None,
        inchi_key: Optional[str] = None,
        ld_id: Optional[str] = None,
        id: Optional[str] = None,
        name: Optional[str] = None,
    ) -> Molecule:
        params = {
            "smiles": smiles,
            "inchi_key": inchi_key,
            "ld_id": ld_id,
            "id": id,
            "name": name,
        }
        if not ld_id:
            params.pop("ld_id")

        molecule = Molecule(**params)

        self.species.append(molecule)

        return molecule

    def add_protein_to_species(
        self,
        sequence: Optional[str] = None,
        organism: Optional[str] = None,
        organism_tax_id: Optional[str] = None,
        ld_id: Optional[str] = None,
        id: Optional[str] = None,
        name: Optional[str] = None,
    ):
        params = {
            "sequence": sequence,
            "organism": organism,
            "organism_tax_id": organism_tax_id,
            "ld_id": ld_id,
            "id": id,
            "name": name,
        }
        if not ld_id:
            params.pop("ld_id")

        protein = Protein(**params)

        self.species.append(protein)

        return protein


class Species(BaseModel):
    model_config: ConfigDict = ConfigDict(  # type: ignore
        validate_assigment=True,
    )  # type: ignore

    id: str
    name: str

    # JSON-LD fields
    ld_id: str = Field(
        serialization_alias="@id", default_factory=lambda: "md:Species/" + str(uuid4())
    )
    ld_type: list[str] = Field(
        serialization_alias="@type",
        default_factory=lambda: [
            "md:Species",
        ],
    )
    ld_context: dict[str, str | dict] = Field(
        serialization_alias="@context",
        default_factory=lambda: {
            "md": "https://github.com/FAIRChemistry/MTPHandler",
        },
    )

    def set_attr_term(
        self,
        attr: str,
        term: str | dict,
        prefix: str | None = None,
        iri: str | None = None,
    ):
        """Sets the term for a given attribute in the JSON-LD object

        Example:
            # Using an IRI term
            >> obj.set_attr_term("name", "http://schema.org/givenName")

            # Using a prefix and term
            >> obj.set_attr_term("name", "schema:givenName", "schema", "http://schema.org")

            # Usinng a dictionary term
            >> obj.set_attr_term("name", {"@id": "http://schema.org/givenName", "@type": "@id"})

        Args:
            attr (str): The attribute to set the term for
            term (str | dict): The term to set for the attribute

        Raises:
            AssertionError: If the attribute is not found in the model
        """

        assert (
            attr in self.model_fields
        ), f"Attribute {attr} not found in {self.__class__.__name__}"

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_context[attr] = term

    def add_type_term(
        self, term: str, prefix: str | None = None, iri: str | None = None
    ):
        """Adds a term to the @type field of the JSON-LD object

        Example:
            # Using a term
            >> obj.add_type_term("https://schema.org/Person")

            # Using a prefixed term
            >> obj.add_type_term("schema:Person", "schema", "https://schema.org/Person")

        Args:
            term (str): The term to add to the @type field
            prefix (str, optional): The prefix to use for the term. Defaults to None.
            iri (str, optional): The IRI to use for the term prefix. Defaults to None.

        Raises:
            ValueError: If prefix is provided but iri is not
            ValueError: If iri is provided but prefix is not
        """

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_type.append(term)


class Molecule(BaseModel):
    model_config: ConfigDict = ConfigDict(  # type: ignore
        validate_assigment=True,
    )  # type: ignore

    smiles: Optional[str] = Field(default=None)
    inchi_key: Optional[str] = Field(default=None)
    id: str
    name: str

    # JSON-LD fields
    ld_id: str = Field(
        serialization_alias="@id", default_factory=lambda: "md:Molecule/" + str(uuid4())
    )
    ld_type: list[str] = Field(
        serialization_alias="@type",
        default_factory=lambda: [
            "md:Molecule",
        ],
    )
    ld_context: dict[str, str | dict] = Field(
        serialization_alias="@context",
        default_factory=lambda: {
            "md": "https://github.com/FAIRChemistry/MTPHandler",
        },
    )

    def set_attr_term(
        self,
        attr: str,
        term: str | dict,
        prefix: str | None = None,
        iri: str | None = None,
    ):
        """Sets the term for a given attribute in the JSON-LD object

        Example:
            # Using an IRI term
            >> obj.set_attr_term("name", "http://schema.org/givenName")

            # Using a prefix and term
            >> obj.set_attr_term("name", "schema:givenName", "schema", "http://schema.org")

            # Usinng a dictionary term
            >> obj.set_attr_term("name", {"@id": "http://schema.org/givenName", "@type": "@id"})

        Args:
            attr (str): The attribute to set the term for
            term (str | dict): The term to set for the attribute

        Raises:
            AssertionError: If the attribute is not found in the model
        """

        assert (
            attr in self.model_fields
        ), f"Attribute {attr} not found in {self.__class__.__name__}"

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_context[attr] = term

    def add_type_term(
        self, term: str, prefix: str | None = None, iri: str | None = None
    ):
        """Adds a term to the @type field of the JSON-LD object

        Example:
            # Using a term
            >> obj.add_type_term("https://schema.org/Person")

            # Using a prefixed term
            >> obj.add_type_term("schema:Person", "schema", "https://schema.org/Person")

        Args:
            term (str): The term to add to the @type field
            prefix (str, optional): The prefix to use for the term. Defaults to None.
            iri (str, optional): The IRI to use for the term prefix. Defaults to None.

        Raises:
            ValueError: If prefix is provided but iri is not
            ValueError: If iri is provided but prefix is not
        """

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_type.append(term)


class Protein(BaseModel):
    model_config: ConfigDict = ConfigDict(  # type: ignore
        validate_assigment=True,
    )  # type: ignore

    sequence: Optional[str] = Field(default=None)
    organism: Optional[str] = Field(default=None)
    organism_tax_id: Optional[int] = Field(default=None)
    id: str
    name: str

    # JSON-LD fields
    ld_id: str = Field(
        serialization_alias="@id", default_factory=lambda: "md:Protein/" + str(uuid4())
    )
    ld_type: list[str] = Field(
        serialization_alias="@type",
        default_factory=lambda: [
            "md:Protein",
        ],
    )
    ld_context: dict[str, str | dict] = Field(
        serialization_alias="@context",
        default_factory=lambda: {
            "md": "https://github.com/FAIRChemistry/MTPHandler",
        },
    )

    def set_attr_term(
        self,
        attr: str,
        term: str | dict,
        prefix: str | None = None,
        iri: str | None = None,
    ):
        """Sets the term for a given attribute in the JSON-LD object

        Example:
            # Using an IRI term
            >> obj.set_attr_term("name", "http://schema.org/givenName")

            # Using a prefix and term
            >> obj.set_attr_term("name", "schema:givenName", "schema", "http://schema.org")

            # Usinng a dictionary term
            >> obj.set_attr_term("name", {"@id": "http://schema.org/givenName", "@type": "@id"})

        Args:
            attr (str): The attribute to set the term for
            term (str | dict): The term to set for the attribute

        Raises:
            AssertionError: If the attribute is not found in the model
        """

        assert (
            attr in self.model_fields
        ), f"Attribute {attr} not found in {self.__class__.__name__}"

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_context[attr] = term

    def add_type_term(
        self, term: str, prefix: str | None = None, iri: str | None = None
    ):
        """Adds a term to the @type field of the JSON-LD object

        Example:
            # Using a term
            >> obj.add_type_term("https://schema.org/Person")

            # Using a prefixed term
            >> obj.add_type_term("schema:Person", "schema", "https://schema.org/Person")

        Args:
            term (str): The term to add to the @type field
            prefix (str, optional): The prefix to use for the term. Defaults to None.
            iri (str, optional): The IRI to use for the term prefix. Defaults to None.

        Raises:
            ValueError: If prefix is provided but iri is not
            ValueError: If iri is provided but prefix is not
        """

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_type.append(term)


class Well(BaseModel):
    model_config: ConfigDict = ConfigDict(  # type: ignore
        validate_assigment=True,
    )  # type: ignore

    id: str
    x_pos: int
    y_pos: int
    ph: Optional[float] = Field(default=None)
    init_conditions: list[InitCondition] = Field(default_factory=list)
    measurements: list[PhotometricMeasurement] = Field(default_factory=list)
    volume: Optional[float] = Field(default=None)
    volume_unit: Optional[UnitDefinition] = Field(default=None)

    # JSON-LD fields
    ld_id: str = Field(
        serialization_alias="@id", default_factory=lambda: "md:Well/" + str(uuid4())
    )
    ld_type: list[str] = Field(
        serialization_alias="@type",
        default_factory=lambda: [
            "md:Well",
        ],
    )
    ld_context: dict[str, str | dict] = Field(
        serialization_alias="@context",
        default_factory=lambda: {
            "md": "https://github.com/FAIRChemistry/MTPHandler",
        },
    )

    def filter_init_conditions(self, **kwargs) -> list[InitCondition]:
        """Filters the init_conditions attribute based on the given kwargs

        Args:
            **kwargs: The attributes to filter by.

        Returns:
            list[InitCondition]: The filtered list of InitCondition objects
        """

        return FilterWrapper[InitCondition](self.init_conditions, **kwargs).filter()

    def filter_measurements(self, **kwargs) -> list[PhotometricMeasurement]:
        """Filters the measurements attribute based on the given kwargs

        Args:
            **kwargs: The attributes to filter by.

        Returns:
            list[PhotometricMeasurement]: The filtered list of PhotometricMeasurement objects
        """

        return FilterWrapper[PhotometricMeasurement](
            self.measurements, **kwargs
        ).filter()

    def set_attr_term(
        self,
        attr: str,
        term: str | dict,
        prefix: str | None = None,
        iri: str | None = None,
    ):
        """Sets the term for a given attribute in the JSON-LD object

        Example:
            # Using an IRI term
            >> obj.set_attr_term("name", "http://schema.org/givenName")

            # Using a prefix and term
            >> obj.set_attr_term("name", "schema:givenName", "schema", "http://schema.org")

            # Usinng a dictionary term
            >> obj.set_attr_term("name", {"@id": "http://schema.org/givenName", "@type": "@id"})

        Args:
            attr (str): The attribute to set the term for
            term (str | dict): The term to set for the attribute

        Raises:
            AssertionError: If the attribute is not found in the model
        """

        assert (
            attr in self.model_fields
        ), f"Attribute {attr} not found in {self.__class__.__name__}"

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_context[attr] = term

    def add_type_term(
        self, term: str, prefix: str | None = None, iri: str | None = None
    ):
        """Adds a term to the @type field of the JSON-LD object

        Example:
            # Using a term
            >> obj.add_type_term("https://schema.org/Person")

            # Using a prefixed term
            >> obj.add_type_term("schema:Person", "schema", "https://schema.org/Person")

        Args:
            term (str): The term to add to the @type field
            prefix (str, optional): The prefix to use for the term. Defaults to None.
            iri (str, optional): The IRI to use for the term prefix. Defaults to None.

        Raises:
            ValueError: If prefix is provided but iri is not
            ValueError: If iri is provided but prefix is not
        """

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_type.append(term)

    def add_to_init_conditions(
        self,
        species_id: str,
        init_conc: float,
        conc_unit: UnitDefinition,
        **kwargs,
    ):
        params = {
            "species_id": species_id,
            "init_conc": init_conc,
            "conc_unit": conc_unit,
        }

        if "id" in kwargs:
            params["id"] = kwargs["id"]

        self.init_conditions.append(InitCondition(**params))

        return self.init_conditions[-1]

    def add_to_measurements(
        self,
        wavelength: float,
        wavelength_unit: UnitDefinition,
        absorption: list[float] = [],
        time: list[float] = [],
        time_unit: Optional[UnitDefinition] = None,
        blank_states: list[BlankState] = [],
        **kwargs,
    ):
        params = {
            "wavelength": wavelength,
            "wavelength_unit": wavelength_unit,
            "absorption": absorption,
            "time": time,
            "time_unit": time_unit,
            "blank_states": blank_states,
        }

        if "id" in kwargs:
            params["id"] = kwargs["id"]

        self.measurements.append(PhotometricMeasurement(**params))

        return self.measurements[-1]

class InitCondition(BaseModel):
    model_config: ConfigDict = ConfigDict(  # type: ignore
        validate_assigment=True,
    )  # type: ignore

    species_id: str
    init_conc: float
    conc_unit: UnitDefinition

    # JSON-LD fields
    ld_id: str = Field(
        serialization_alias="@id",
        default_factory=lambda: "md:PhotometricMeasurement/" + str(uuid4()),
    )
    ld_type: list[str] = Field(
        serialization_alias="@type",
        default_factory=lambda: [
            "md:InitCondition",
        ],
    )
    ld_context: dict[str, str | dict] = Field(
        serialization_alias="@context",
        default_factory=lambda: {
            "md": "https://github.com/FAIRChemistry/MTPHandler",
        },
    )


class PhotometricMeasurement(BaseModel):
    model_config: ConfigDict = ConfigDict(  # type: ignore
        validate_assigment=True,
    )  # type: ignore

    wavelength: float
    wavelength_unit: UnitDefinition
    absorption: list[float] = Field(default_factory=list)
    time: list[float] = Field(default_factory=list)
    time_unit: Optional[UnitDefinition] = Field(default=None)
    blank_states: list[BlankState] = Field(default_factory=list)

    # JSON-LD fields
    ld_id: str = Field(
        serialization_alias="@id",
        default_factory=lambda: "md:PhotometricMeasurement/" + str(uuid4()),
    )
    ld_type: list[str] = Field(
        serialization_alias="@type",
        default_factory=lambda: [
            "md:PhotometricMeasurement",
        ],
    )
    ld_context: dict[str, str | dict] = Field(
        serialization_alias="@context",
        default_factory=lambda: {
            "md": "https://github.com/FAIRChemistry/MTPHandler",
        },
    )

    def filter_blank_states(self, **kwargs) -> list[BlankState]:
        """Filters the blank_states attribute based on the given kwargs

        Args:
            **kwargs: The attributes to filter by.

        Returns:
            list[BlankState]: The filtered list of BlankState objects
        """

        return FilterWrapper[BlankState](self.blank_states, **kwargs).filter()

    def set_attr_term(
        self,
        attr: str,
        term: str | dict,
        prefix: str | None = None,
        iri: str | None = None,
    ):
        """Sets the term for a given attribute in the JSON-LD object

        Example:
            # Using an IRI term
            >> obj.set_attr_term("name", "http://schema.org/givenName")

            # Using a prefix and term
            >> obj.set_attr_term("name", "schema:givenName", "schema", "http://schema.org")

            # Usinng a dictionary term
            >> obj.set_attr_term("name", {"@id": "http://schema.org/givenName", "@type": "@id"})

        Args:
            attr (str): The attribute to set the term for
            term (str | dict): The term to set for the attribute

        Raises:
            AssertionError: If the attribute is not found in the model
        """

        assert (
            attr in self.model_fields
        ), f"Attribute {attr} not found in {self.__class__.__name__}"

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_context[attr] = term

    def add_type_term(
        self, term: str, prefix: str | None = None, iri: str | None = None
    ):
        """Adds a term to the @type field of the JSON-LD object

        Example:
            # Using a term
            >> obj.add_type_term("https://schema.org/Person")

            # Using a prefixed term
            >> obj.add_type_term("schema:Person", "schema", "https://schema.org/Person")

        Args:
            term (str): The term to add to the @type field
            prefix (str, optional): The prefix to use for the term. Defaults to None.
            iri (str, optional): The IRI to use for the term prefix. Defaults to None.

        Raises:
            ValueError: If prefix is provided but iri is not
            ValueError: If iri is provided but prefix is not
        """

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_type.append(term)

    def add_to_blank_states(
        self,
        species_id: str,
        contributes_to_signal: bool = True,
        **kwargs,
    ):
        params = {
            "species_id": species_id,
            "contributes_to_signal": contributes_to_signal,
        }

        if "id" in kwargs:
            params["id"] = kwargs["id"]

        self.blank_states.append(BlankState(**params))

        return self.blank_states[-1]


class spec(BaseModel):
    model_config: ConfigDict = ConfigDict(  # type: ignore
        validate_assigment=True,
    )  # type: ignore

    species_id: str
    init_conc: float
    conc_unit: UnitDefinition

    # JSON-LD fields
    ld_id: str = Field(
        serialization_alias="@id",
        default_factory=lambda: "md:InitCondition/" + str(uuid4()),
    )
    ld_type: list[str] = Field(
        serialization_alias="@type",
        default_factory=lambda: [
            "md:InitCondition",
        ],
    )
    ld_context: dict[str, str | dict] = Field(
        serialization_alias="@context",
        default_factory=lambda: {
            "md": "https://github.com/FAIRChemistry/MTPHandler",
        },
    )

    def set_attr_term(
        self,
        attr: str,
        term: str | dict,
        prefix: str | None = None,
        iri: str | None = None,
    ):
        """Sets the term for a given attribute in the JSON-LD object

        Example:
            # Using an IRI term
            >> obj.set_attr_term("name", "http://schema.org/givenName")

            # Using a prefix and term
            >> obj.set_attr_term("name", "schema:givenName", "schema", "http://schema.org")

            # Usinng a dictionary term
            >> obj.set_attr_term("name", {"@id": "http://schema.org/givenName", "@type": "@id"})

        Args:
            attr (str): The attribute to set the term for
            term (str | dict): The term to set for the attribute

        Raises:
            AssertionError: If the attribute is not found in the model
        """

        assert (
            attr in self.model_fields
        ), f"Attribute {attr} not found in {self.__class__.__name__}"

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_context[attr] = term

    def add_type_term(
        self, term: str, prefix: str | None = None, iri: str | None = None
    ):
        """Adds a term to the @type field of the JSON-LD object

        Example:
            # Using a term
            >> obj.add_type_term("https://schema.org/Person")

            # Using a prefixed term
            >> obj.add_type_term("schema:Person", "schema", "https://schema.org/Person")

        Args:
            term (str): The term to add to the @type field
            prefix (str, optional): The prefix to use for the term. Defaults to None.
            iri (str, optional): The IRI to use for the term prefix. Defaults to None.

        Raises:
            ValueError: If prefix is provided but iri is not
            ValueError: If iri is provided but prefix is not
        """

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_type.append(term)


class BlankState(BaseModel):
    model_config: ConfigDict = ConfigDict(  # type: ignore
        validate_assigment=True,
    )  # type: ignore

    species_id: str
    contributes_to_signal: bool = True

    # JSON-LD fields
    ld_id: str = Field(
        serialization_alias="@id",
        default_factory=lambda: "md:BlankState/" + str(uuid4()),
    )
    ld_type: list[str] = Field(
        serialization_alias="@type",
        default_factory=lambda: [
            "md:BlankState",
        ],
    )
    ld_context: dict[str, str | dict] = Field(
        serialization_alias="@context",
        default_factory=lambda: {
            "md": "https://github.com/FAIRChemistry/MTPHandler",
        },
    )

    def set_attr_term(
        self,
        attr: str,
        term: str | dict,
        prefix: str | None = None,
        iri: str | None = None,
    ):
        """Sets the term for a given attribute in the JSON-LD object

        Example:
            # Using an IRI term
            >> obj.set_attr_term("name", "http://schema.org/givenName")

            # Using a prefix and term
            >> obj.set_attr_term("name", "schema:givenName", "schema", "http://schema.org")

            # Usinng a dictionary term
            >> obj.set_attr_term("name", {"@id": "http://schema.org/givenName", "@type": "@id"})

        Args:
            attr (str): The attribute to set the term for
            term (str | dict): The term to set for the attribute

        Raises:
            AssertionError: If the attribute is not found in the model
        """

        assert (
            attr in self.model_fields
        ), f"Attribute {attr} not found in {self.__class__.__name__}"

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_context[attr] = term

    def add_type_term(
        self, term: str, prefix: str | None = None, iri: str | None = None
    ):
        """Adds a term to the @type field of the JSON-LD object

        Example:
            # Using a term
            >> obj.add_type_term("https://schema.org/Person")

            # Using a prefixed term
            >> obj.add_type_term("schema:Person", "schema", "https://schema.org/Person")

        Args:
            term (str): The term to add to the @type field
            prefix (str, optional): The prefix to use for the term. Defaults to None.
            iri (str, optional): The IRI to use for the term prefix. Defaults to None.

        Raises:
            ValueError: If prefix is provided but iri is not
            ValueError: If iri is provided but prefix is not
        """

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_type.append(term)


class UnitDefinition(BaseModel):
    model_config: ConfigDict = ConfigDict(  # type: ignore
        validate_assigment=True, use_enum_values=True
    )  # type: ignore

    id: Optional[str] = Field(default=None)
    name: Optional[str] = Field(default=None)
    base_units: list[BaseUnit] = Field(default_factory=list)

    # JSON-LD fields
    ld_id: str = Field(
        serialization_alias="@id",
        default_factory=lambda: "md:UnitDefinition/" + str(uuid4()),
    )
    ld_type: list[str] = Field(
        serialization_alias="@type",
        default_factory=lambda: [
            "md:UnitDefinition",
        ],
    )
    ld_context: dict[str, str | dict] = Field(
        serialization_alias="@context",
        default_factory=lambda: {
            "md": "https://github.com/FAIRChemistry/MTPHandler",
        },
    )

    def filter_base_units(self, **kwargs) -> list[BaseUnit]:
        """Filters the base_units attribute based on the given kwargs

        Args:
            **kwargs: The attributes to filter by.

        Returns:
            list[BaseUnit]: The filtered list of BaseUnit objects
        """

        return FilterWrapper[BaseUnit](self.base_units, **kwargs).filter()

    def set_attr_term(
        self,
        attr: str,
        term: str | dict,
        prefix: str | None = None,
        iri: str | None = None,
    ):
        """Sets the term for a given attribute in the JSON-LD object

        Example:
            # Using an IRI term
            >> obj.set_attr_term("name", "http://schema.org/givenName")

            # Using a prefix and term
            >> obj.set_attr_term("name", "schema:givenName", "schema", "http://schema.org")

            # Usinng a dictionary term
            >> obj.set_attr_term("name", {"@id": "http://schema.org/givenName", "@type": "@id"})

        Args:
            attr (str): The attribute to set the term for
            term (str | dict): The term to set for the attribute

        Raises:
            AssertionError: If the attribute is not found in the model
        """

        assert (
            attr in self.model_fields
        ), f"Attribute {attr} not found in {self.__class__.__name__}"

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_context[attr] = term

    def add_type_term(
        self, term: str, prefix: str | None = None, iri: str | None = None
    ):
        """Adds a term to the @type field of the JSON-LD object

        Example:
            # Using a term
            >> obj.add_type_term("https://schema.org/Person")

            # Using a prefixed term
            >> obj.add_type_term("schema:Person", "schema", "https://schema.org/Person")

        Args:
            term (str): The term to add to the @type field
            prefix (str, optional): The prefix to use for the term. Defaults to None.
            iri (str, optional): The IRI to use for the term prefix. Defaults to None.

        Raises:
            ValueError: If prefix is provided but iri is not
            ValueError: If iri is provided but prefix is not
        """

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_type.append(term)

    def add_to_base_units(
        self,
        kind: UnitType,
        exponent: int,
        multiplier: Optional[float] = None,
        scale: Optional[float] = None,
        **kwargs,
    ):
        params = {
            "kind": kind,
            "exponent": exponent,
            "multiplier": multiplier,
            "scale": scale,
        }

        if "id" in kwargs:
            params["id"] = kwargs["id"]

        self.base_units.append(BaseUnit(**params))

        return self.base_units[-1]


class BaseUnit(BaseModel):
    model_config: ConfigDict = ConfigDict(  # type: ignore
        validate_assigment=True, use_enum_values=True
    )  # type: ignore

    kind: UnitType
    exponent: int
    multiplier: Optional[float] = Field(default=None)
    scale: Optional[float] = Field(default=None)

    # JSON-LD fields
    ld_id: str = Field(
        serialization_alias="@id", default_factory=lambda: "md:BaseUnit/" + str(uuid4())
    )
    ld_type: list[str] = Field(
        serialization_alias="@type",
        default_factory=lambda: [
            "md:BaseUnit",
        ],
    )
    ld_context: dict[str, str | dict] = Field(
        serialization_alias="@context",
        default_factory=lambda: {
            "md": "https://github.com/FAIRChemistry/MTPHandler",
        },
    )

    def set_attr_term(
        self,
        attr: str,
        term: str | dict,
        prefix: str | None = None,
        iri: str | None = None,
    ):
        """Sets the term for a given attribute in the JSON-LD object

        Example:
            # Using an IRI term
            >> obj.set_attr_term("name", "http://schema.org/givenName")

            # Using a prefix and term
            >> obj.set_attr_term("name", "schema:givenName", "schema", "http://schema.org")

            # Usinng a dictionary term
            >> obj.set_attr_term("name", {"@id": "http://schema.org/givenName", "@type": "@id"})

        Args:
            attr (str): The attribute to set the term for
            term (str | dict): The term to set for the attribute

        Raises:
            AssertionError: If the attribute is not found in the model
        """

        assert (
            attr in self.model_fields
        ), f"Attribute {attr} not found in {self.__class__.__name__}"

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_context[attr] = term

    def add_type_term(
        self, term: str, prefix: str | None = None, iri: str | None = None
    ):
        """Adds a term to the @type field of the JSON-LD object

        Example:
            # Using a term
            >> obj.add_type_term("https://schema.org/Person")

            # Using a prefixed term
            >> obj.add_type_term("schema:Person", "schema", "https://schema.org/Person")

        Args:
            term (str): The term to add to the @type field
            prefix (str, optional): The prefix to use for the term. Defaults to None.
            iri (str, optional): The IRI to use for the term prefix. Defaults to None.

        Raises:
            ValueError: If prefix is provided but iri is not
            ValueError: If iri is provided but prefix is not
        """

        if prefix:
            validate_prefix(term, prefix)

        add_namespace(self, prefix, iri)
        self.ld_type.append(term)


class UnitType(Enum):
    AMPERE = "ampere"
    AVOGADRO = "avogadro"
    BECQUEREL = "becquerel"
    CANDELA = "candela"
    CELSIUS = "celsius"
    COULOMB = "coulomb"
    DIMENSIONLESS = "dimensionless"
    FARAD = "farad"
    GRAM = "gram"
    GRAY = "gray"
    HENRY = "henry"
    HERTZ = "hertz"
    ITEM = "item"
    JOULE = "joule"
    KATAL = "katal"
    KELVIN = "kelvin"
    KILOGRAM = "kilogram"
    LITRE = "litre"
    LUMEN = "lumen"
    LUX = "lux"
    METRE = "metre"
    MOLE = "mole"
    NEWTON = "newton"
    OHM = "ohm"
    PASCAL = "pascal"
    RADIAN = "radian"
    SECOND = "second"
    SIEMENS = "siemens"
    SIEVERT = "sievert"
    STERADIAN = "steradian"
    TESLA = "tesla"
    VOLT = "volt"
    WATT = "watt"
    WEBER = "weber"
