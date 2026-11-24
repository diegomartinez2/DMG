import time
import cv2
import numpy as np

cap = cv2.VideoCapture(0)

last_update = 0.0
hold_duration = 0.0
display_frame = None

while True:
    now = time.monotonic()

    # Decidir si toca capturar y generar un nuevo frame
    if display_frame is None or (now - last_update) >= hold_duration:
        ret, frame = cap.read()
        if not ret:
            break

        h, w, _ = frame.shape
        noisy = frame.astype(np.float32)

        # 1. Ruido granular (estática)
        noise = np.random.normal(0, 28, (h, w, 1)).astype(np.float32)
        # Strong:
        # noise = np.random.normal(0, 60, (h, w, 1)).astype(np.float32)
        # Weak:
        # noise = np.random.normal(0, 10, (h, w, 1)).astype(np.float32)
        noisy = np.clip(noisy + noise, 0, 255).astype(np.uint8)

        # 2. Líneas horizontales de interferencia
        for _ in range(np.random.randint(2, 7)):
            y = np.random.randint(0, h - 1)
            band_h = np.random.randint(1, 5)
            shift = np.random.randint(-25, 25)
            noisy[y:y + band_h] = np.roll(noisy[y:y + band_h], shift, axis=1)

        # Strong:
        '''
        for _ in range(np.random.randint(10, 22)):
            y = np.random.randint(0, h - 1)
            band_h = np.random.randint(3, 12)
            shift = np.random.randint(-70, 70)
            noisy[y:y + band_h] = np.roll(noisy[y:y + band_h], shift, axis=1)
        '''
        # Weak:
        '''
        for _ in range(np.random.randint(1, 3)):
            y = np.random.randint(0, h - 1)
            band_h = np.random.randint(1, 3)
            shift = np.random.randint(-8, 8)
            noisy[y:y + band_h] = np.roll(noisy[y:y + band_h], shift, axis=1)
        '''

        display_frame = noisy
        last_update = now
        # Duración aleatoria del frame entre 0 y 3 segundos
        hold_duration = np.random.uniform(0.0, 3.0)

    if display_frame is not None:
        cv2.imshow("Webcam con interferencia", display_frame)

    if cv2.waitKey(15) & 0xFF == ord("q"):
        break

cap.release()
cv2.destroyAllWindows()
