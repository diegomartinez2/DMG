#!/usr/bin/env python
import asyncio
from pyhanko.keys import SimpleSigner
from pyhanko.pdf_utils.incremental_writer import IncrementalPdfFileWriter
from pyhanko.sign import fields, PdfSigner
from pyhanko.sign.timestamps import TimeStamper
from pyhanko.sign.fields import SigFieldSpec
from pyhanko.sign.ades.api import CAdESSignedData
from pyhanko.sign.diff_analysis import summarize_and_validate_pdf_with_diff

# 1. Rutas a tus archivos
PDF_INPUT = 'documento.pdf'
PDF_OUTPUT = 'documento_firmado_visible.pdf'
CERT_FILE = 'cert.pem'#'mi_certificado.pem'  # o .crt
KEY_FILE = 'myKey.pem'#'mi_clave_privada.pem' # o .key
KEY_PASS = b'tu_contrasena' # Si tu clave tiene contraseña
CA_CHAIN_FILES = ('tu_ca_intermedia.pem', 'tu_ca_raiz.pem') # Opcional: cadena de CA

async def sign_pdf_visible():
    # 2. Cargar el firmante
    signer = SimpleSigner.load(
        KEY_FILE,
        CERT_FILE,
        ca_chain_files=CA_CHAIN_FILES,
        key_passphrase=KEY_PASS
    )

    with open(PDF_INPUT, 'rb') as inf:
        w = IncrementalPdfFileWriter(inf)

        # 3. Crear una especificación de campo de firma visible
        # box=(x1, y1, x2, y2)
        # x1, y1: esquina inferior izquierda
        # x2, y2: esquina superior derecha
        # on_page: índice de la página (0 para la primera, -1 para la última)
        sig_field_spec = SigFieldSpec(
            'MiFirmaVisible',
            box=(50, 50, 250, 150),  # Ajusta estas coordenadas y tamaño
            on_page=0 # Página 1 (índice 0)
        )

        # Añadir el campo de firma (esto también se puede hacer por separado)
        fields.append_signature_field(w, sig_field_spec)

        # 4. Personalizar la apariencia del sello (opcional)
        from pyhanko.stamp import TextStampStyle, QRStampStyle
        from pyhanko.pdf_utils import text, images
        from pyhanko.fonts import opentype

        # Ejemplo de sello de texto:
        # Aquí puedes usar %(signer)s, %(ts)s, %(reason)s, %(location)s, etc.
        # pyHanko los interpolará automáticamente
        text_stamp_style = TextStampStyle(
            stamp_text='Firmado por: %(signer)s\nFecha: %(ts)s',
            text_box_style=text.TextBoxStyle(
                font=opentype.GlyphAccumulatorFactory('ruta/a/tu_fuente.ttf') # Opcional: usar una fuente personalizada
            ),
            # background=images.PdfImage.from_file('ruta/a/imagen_fondo.png') # Opcional: añadir una imagen de fondo
        )

        # Ejemplo de sello con código QR (útil para URLs de verificación)
        # qr_stamp_style = QRStampStyle(
        #     stamp_text='Verificar: %(url)s\nFirmado por: %(signer)s',
        #     url_field='url_de_verificacion' # Necesitarías añadir esto a tus metadatos
        # )

        # 5. Configurar los metadatos de la firma
        meta = PdfSigner.recursively_get_sig_metadata(
            w, field_name='MiFirmaVisible',
            # Opcional: Puedes añadir más metadatos para que aparezcan en el sello
            # reason='Aprobación del documento',
            # location='Cagliari, Italia',
            # contact_info='info@ejemplo.com'
        )

        # Inicializar el firmador con el estilo de sello deseado
        pdf_signer = PdfSigner(
            meta,
            signer=signer,
            stamp_style=text_stamp_style # Usa el estilo de sello que definiste
        )

        # Firmar el documento
        output_buffer = IncrementalPdfFileWriter(inf)
        # Asegúrate de usar 'w' (IncrementalPdfFileWriter) que contiene el campo de firma
        async with open(PDF_OUTPUT, 'wb') as outf:
            await pdf_signer.sign_pdf(
                w, output_buffer, output=outf,
                # timestamper=TimeStamper(url='http://timestamp.digicert.com') # Opcional: servidor de sello de tiempo
            )

    print(f"Documento firmado con firma visible guardado en: {PDF_OUTPUT}")

# Ejecutar la función asíncrona
if __name__ == '__main__':
    # Asegúrate de que los archivos de certificado y clave existan y sean válidos.
    # Si usas una fuente personalizada, asegúrate de que la ruta sea correcta.
    asyncio.run(sign_pdf_visible())
