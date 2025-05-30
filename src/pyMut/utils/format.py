def formatear_rs(cadena):
    # Separamos la cadena por '|' para obtener cada código
    codigos = cadena.split('|')
    # Eliminamos el prefijo "rs" de cada código
    codigos_solo_numeros = [codigo[2:] if codigo.startswith("rs") else codigo for codigo in codigos]

    return '|'.join(codigos_solo_numeros)

def formatear_chr(string:str):
    """
    Formatea el cromosoma a un formato estándar.
    Si aparece como 1, lo convierte a chr1.
    Convierte "23" y "24" a "X" e "Y" respectivamente.
    """
    if string == "23":
        return "X"
    elif string == "24":
        return "Y"
    elif string.startswith("chr"):
        return string
    else:
        return "chr" + string